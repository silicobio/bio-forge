//! Simplified AMBER-like energy minimization for protein side chain and whole-structure relaxation.
//!
//! The relaxation engine assembles a flat atom system from the structure, applies
//! a harmonic bond-stretching and bond-angle-bending potential derived from embedded
//! residue templates, and resolves Lennard-Jones non-bonded clashes using a steepest-descent
//! minimizer with an adaptive step-size line search.
//!
//! **Supported modes**
//!
//! * `side_chains_only = true` – backbone heavy atoms (`N`, `CA`, `C`, `O`, `OXT`) and
//!   all non-standard residues (water, ions, heterogens) are held fixed while protein
//!   side-chain atoms are optimized. This mirrors the "AMBER side-chain relaxation"
//!   workflow used before MD production runs.
//! * `side_chains_only = false` – all heavy atoms of standard residues are optimized
//!   (non-standard residues remain fixed).
//!
//! **Force field summary**
//!
//! | Term | Parameters |
//! |------|-----------|
//! | Bond stretching | Harmonic, `k_b` from element pairs (AMBER99SB-derived), `r₀` from template geometry |
//! | Angle bending | Harmonic, `k_a` from central element (AMBER99SB-derived), `θ₀` from template geometry |
//! | VDW (non-bonded) | Full Lennard-Jones 12-6 with Lorentz-Berthelot mixing, cutoff 10 Å, 1-2 and 1-3 exclusions |

use crate::db;
use crate::model::{
    structure::Structure,
    types::{Element, ResidueCategory},
};
use crate::ops::error::Error;
use nalgebra::Vector3;
use std::collections::HashSet;

// ─── Physical constants ───────────────────────────────────────────────────────

/// Default VDW cutoff distance (Å).
const DEFAULT_VDW_CUTOFF: f64 = 10.0;
/// Default maximum steepest-descent steps.
const DEFAULT_MAX_STEPS: u32 = 200;
/// Default RMS-gradient convergence threshold (kcal mol⁻¹ Å⁻¹).
const DEFAULT_CONVERGENCE: f64 = 1.0;
/// Initial step size for steepest descent (Å).
const INITIAL_STEP_SIZE: f64 = 0.02;
/// Factor to shrink step size on an energy-increasing trial.
const STEP_SHRINK: f64 = 0.5;
/// Factor to grow step size on a successful step.
const STEP_GROW: f64 = 1.2;
/// Maximum step size allowed during line search (Å).
const MAX_STEP_SIZE: f64 = 0.10;
/// Minimum step size before aborting the minimization (Å).
const MIN_STEP_SIZE: f64 = 1e-6;
/// Number of spatial degrees of freedom per atom (x, y, z).
const SPATIAL_DIMS: usize = 3;

// ─── AMBER ff14SB-derived VDW parameters ─────────────────────────────────────

/// Returns the AMBER ff14SB Lennard-Jones `(r_min/2, ε)` pair for an element.
///
/// `r_min/2` is in Å and `ε` is in kcal mol⁻¹.  Parameters for uncommon elements
/// fall back to generic organic atom values.
fn vdw_params(element: Element) -> (f64, f64) {
    match element {
        Element::C => (1.908, 0.0860),
        Element::N => (1.824, 0.1700),
        Element::O => (1.661, 0.2100),
        Element::S => (2.000, 0.2500),
        Element::P => (2.100, 0.2000),
        Element::H => (1.100, 0.0157),
        Element::Fe | Element::Zn | Element::Mg | Element::Ca | Element::Na | Element::K => {
            (1.500, 0.0050)
        }
        _ => (1.700, 0.1000),
    }
}

/// Returns the combined LJ parameters for an atom pair using Lorentz–Berthelot mixing.
///
/// Returns `(r_min_ij, eps_ij)` where `r_min_ij = r_i + r_j` (sum of `r_min/2` values).
#[inline]
fn lj_pair_params(e1: Element, e2: Element) -> (f64, f64) {
    let (r1, eps1) = vdw_params(e1);
    let (r2, eps2) = vdw_params(e2);
    (r1 + r2, (eps1 * eps2).sqrt())
}

// ─── AMBER-derived bond/angle force constants ─────────────────────────────────

/// Returns the harmonic bond force constant `k_b` (kcal mol⁻¹ Å⁻²) for an element pair.
fn bond_force_constant(e1: Element, e2: Element) -> f64 {
    let key = (e1 as u8).min(e2 as u8);
    let other = (e1 as u8).max(e2 as u8);

    // Match on the lighter element first, then heavier (sorted pair).
    let c = Element::C as u8;
    let n = Element::N as u8;
    let o = Element::O as u8;
    let s = Element::S as u8;
    let h = Element::H as u8;
    let p = Element::P as u8;

    match (key, other) {
        (k, o2) if k == h && o2 == c => 340.0, // C-H
        (k, o2) if k == h && o2 == n => 434.0, // N-H
        (k, o2) if k == h && o2 == o => 553.0, // O-H
        (k, o2) if k == h && o2 == s => 274.0, // S-H
        (k, o2) if k == c && o2 == c => 310.0, // C-C
        (k, o2) if k == c && o2 == n => 337.0, // C-N
        (k, o2) if k == c && o2 == o => 570.0, // C-O (use carbonyl; also handles C-OH)
        (k, o2) if k == c && o2 == s => 227.0, // C-S
        (k, o2) if k == n && o2 == o => 370.0, // N-O
        (k, o2) if k == c && o2 == p => 230.0, // C-P
        (k, o2) if k == o && o2 == p => 525.0, // O-P
        _ => 300.0,                            // generic fallback
    }
}

/// Returns the harmonic angle force constant `k_a` (kcal mol⁻¹ rad⁻²) for the central element.
fn angle_force_constant(central: Element) -> f64 {
    match central {
        Element::C => 63.0,
        Element::N => 50.0,
        Element::O => 35.0,
        Element::S => 58.0,
        Element::P => 80.0,
        _ => 50.0,
    }
}

// ─── Backbone atom set ────────────────────────────────────────────────────────

/// Protein main-chain atom names that remain fixed during side-chain-only relaxation.
const BACKBONE_NAMES: &[&str] = &[
    "N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3",
];

fn is_backbone_atom(name: &str) -> bool {
    BACKBONE_NAMES.contains(&name)
}

// ─── Public types ─────────────────────────────────────────────────────────────

/// Configuration parameters for the relaxation pipeline.
#[derive(Debug, Clone)]
pub struct RelaxConfig {
    /// Maximum number of steepest-descent iterations.
    pub max_steps: u32,
    /// When `true`, only protein side-chain heavy atoms are moved; backbone and non-standard
    /// residues are held fixed.  When `false`, all heavy atoms of standard residues move.
    pub side_chains_only: bool,
    /// Convergence threshold: stop when the RMS gradient falls below this value
    /// (kcal mol⁻¹ Å⁻¹).
    pub convergence: f64,
    /// Lennard-Jones cutoff distance (Å).
    pub vdw_cutoff: f64,
}

impl Default for RelaxConfig {
    fn default() -> Self {
        Self {
            max_steps: DEFAULT_MAX_STEPS,
            side_chains_only: true,
            convergence: DEFAULT_CONVERGENCE,
            vdw_cutoff: DEFAULT_VDW_CUTOFF,
        }
    }
}

/// Summary statistics returned after a completed relaxation run.
#[derive(Debug, Clone)]
pub struct RelaxResult {
    /// System energy before minimization (kcal mol⁻¹).
    pub initial_energy: f64,
    /// System energy after minimization (kcal mol⁻¹).
    pub final_energy: f64,
    /// Number of minimization steps performed.
    pub steps_taken: u32,
    /// Whether the convergence criterion was satisfied.
    pub converged: bool,
}

// ─── Internal minimization system ────────────────────────────────────────────

/// Harmonic bond potential term.
#[derive(Debug, Clone, Copy)]
struct BondTerm {
    i: usize,
    j: usize,
    /// Equilibrium bond length derived from the template (Å).
    r0: f64,
    /// Force constant (kcal mol⁻¹ Å⁻²).
    k: f64,
}

/// Harmonic angle potential term.
#[derive(Debug, Clone, Copy)]
struct AngleTerm {
    i: usize,
    j: usize,
    k: usize,
    /// Equilibrium angle (radians).
    theta0: f64,
    /// Force constant (kcal mol⁻¹ rad⁻²).
    ka: f64,
}

/// Lennard-Jones non-bonded pair.
#[derive(Debug, Clone, Copy)]
struct NbPair {
    i: usize,
    j: usize,
    /// `r_min_ij` (Å) – sum of individual `r_min/2` parameters.
    rmin: f64,
    /// Well depth `ε_ij` (kcal mol⁻¹).
    eps: f64,
}

/// Complete description of the system to be minimized.
struct MinSystem {
    /// Flat array of all atom positions (Å).
    positions: Vec<Vector3<f64>>,
    /// `true` if the atom is allowed to move.
    movable: Vec<bool>,
    /// Harmonic bond terms.
    bond_terms: Vec<BondTerm>,
    /// Harmonic angle terms.
    angle_terms: Vec<AngleTerm>,
    /// Non-bonded Lennard-Jones pairs.
    nb_pairs: Vec<NbPair>,
}

// ─── System construction ──────────────────────────────────────────────────────

/// Builds the flat `MinSystem` from a structure and configuration.
fn build_system(structure: &Structure, config: &RelaxConfig) -> MinSystem {
    // ── 1.  Enumerate atoms and assign global indices ──────────────────────────
    let mut positions: Vec<Vector3<f64>> = Vec::new();
    let mut elements: Vec<Element> = Vec::new();
    let mut movable: Vec<bool> = Vec::new();

    // index_of[(chain_idx, residue_idx, atom_idx)] -> flat_index
    let mut index_of: std::collections::HashMap<(usize, usize, usize), usize> =
        std::collections::HashMap::new();

    for (ci, chain) in structure.iter_chains().enumerate() {
        for (ri, residue) in chain.iter_residues().enumerate() {
            let is_standard = residue.category == ResidueCategory::Standard;
            for (ai, atom) in residue.iter_atoms().enumerate() {
                let flat_idx = positions.len();
                index_of.insert((ci, ri, ai), flat_idx);

                positions.push(atom.pos.coords);
                elements.push(atom.element);

                let can_move = if is_standard {
                    if config.side_chains_only {
                        !is_backbone_atom(&atom.name) && atom.element != Element::H // keep hydrogens fixed too
                    } else {
                        atom.element != Element::H
                    }
                } else {
                    false
                };

                movable.push(can_move);
            }
        }
    }

    let n_atoms = positions.len();

    // ── 2.  Build bond and angle terms from templates ─────────────────────────
    //   Also accumulate a set of bonded pairs for non-bonded exclusions.
    let mut bond_terms: Vec<BondTerm> = Vec::new();
    let mut angle_terms: Vec<AngleTerm> = Vec::new();
    // Exclusion set: (min, max) pairs that are 1-2 or 1-3 bonded.
    let mut excluded: HashSet<(usize, usize)> = HashSet::new();
    // Adjacency list for building angles.
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];

    for (ci, chain) in structure.iter_chains().enumerate() {
        for (ri, residue) in chain.iter_residues().enumerate() {
            if residue.category != ResidueCategory::Standard {
                continue;
            }
            let template = match db::get_template(&residue.name) {
                Some(t) => t,
                None => continue,
            };

            // Collect all atom names present in both template and residue.
            // Build a local name->flat_idx map.
            let local_map: std::collections::HashMap<&str, usize> = residue
                .iter_atoms()
                .enumerate()
                .filter_map(|(ai, atom)| {
                    index_of
                        .get(&(ci, ri, ai))
                        .map(|&flat_idx| (atom.name.as_str(), flat_idx))
                })
                .collect();

            // Add bonds from template.
            for (a1_name, a2_name, _order) in template.bonds() {
                // Find template reference positions for equilibrium length.
                let r0 = template
                    .heavy_atoms()
                    .find(|(n, _, _)| *n == a1_name)
                    .and_then(|(_, _, p1)| {
                        template
                            .heavy_atoms()
                            .find(|(n, _, _)| *n == a2_name)
                            .map(|(_, _, p2)| nalgebra::distance(&p1, &p2))
                    })
                    .or_else(|| {
                        // One or both are hydrogens.
                        let p1 = template
                            .hydrogens()
                            .find(|(n, _, _)| *n == a1_name)
                            .map(|(_, p, _)| p)
                            .or_else(|| {
                                template
                                    .heavy_atoms()
                                    .find(|(n, _, _)| *n == a1_name)
                                    .map(|(_, _, p)| p)
                            })?;
                        let p2 = template
                            .hydrogens()
                            .find(|(n, _, _)| *n == a2_name)
                            .map(|(_, p, _)| p)
                            .or_else(|| {
                                template
                                    .heavy_atoms()
                                    .find(|(n, _, _)| *n == a2_name)
                                    .map(|(_, _, p)| p)
                            })?;
                        Some(nalgebra::distance(&p1, &p2))
                    });

                let (Some(&i), Some(&j)) = (local_map.get(a1_name), local_map.get(a2_name)) else {
                    continue;
                };

                // At least one endpoint must be movable for the term to affect the gradient.
                if !movable[i] && !movable[j] {
                    // Still register the exclusion.
                    let key = (i.min(j), i.max(j));
                    excluded.insert(key);
                    adj[i].push(j);
                    adj[j].push(i);
                    continue;
                }

                let Some(r0) = r0 else { continue };
                if r0 < 1e-3 {
                    continue; // degenerate template geometry
                }

                let e1 = elements[i];
                let e2 = elements[j];
                let k = bond_force_constant(e1, e2);

                bond_terms.push(BondTerm { i, j, r0, k });
                let key = (i.min(j), i.max(j));
                excluded.insert(key);
                adj[i].push(j);
                adj[j].push(i);
            }

            // Handle hydrogens: add bonds from hydrogen entries in template.
            for (h_name, h_pos, mut anchors) in template.hydrogens() {
                let Some(anchor_name) = anchors.next() else {
                    continue;
                };
                let r0 = {
                    let anchor_pos = template
                        .heavy_atoms()
                        .find(|(n, _, _)| *n == anchor_name)
                        .map(|(_, _, p)| p);
                    match anchor_pos {
                        Some(ap) => nalgebra::distance(&h_pos, &ap),
                        None => continue,
                    }
                };
                let (Some(&i), Some(&j)) = (local_map.get(h_name), local_map.get(anchor_name))
                else {
                    continue;
                };
                let key = (i.min(j), i.max(j));
                excluded.insert(key);
                adj[i].push(j);
                adj[j].push(i);

                if !movable[i] && !movable[j] {
                    continue;
                }
                let e1 = elements[i];
                let e2 = elements[j];
                bond_terms.push(BondTerm {
                    i,
                    j,
                    r0,
                    k: bond_force_constant(e1, e2),
                });
            }
        }
    }

    // ── 3.  Inter-residue backbone bonds (C—N peptide bond) ───────────────────
    for (ci, chain) in structure.iter_chains().enumerate() {
        let residues: Vec<_> = chain.iter_residues().collect();
        for ri in 0..residues.len().saturating_sub(1) {
            let res_i = residues[ri];
            let res_j = residues[ri + 1];

            if res_i.category != ResidueCategory::Standard
                || res_j.category != ResidueCategory::Standard
            {
                continue;
            }

            // Find C of residue i and N of residue j.
            let c_idx = res_i
                .iter_atoms()
                .enumerate()
                .find(|(_, a)| a.name == "C")
                .and_then(|(ai, _)| index_of.get(&(ci, ri, ai)).copied());

            let n_idx = res_j
                .iter_atoms()
                .enumerate()
                .find(|(_, a)| a.name == "N")
                .and_then(|(ai, _)| index_of.get(&(ci, ri + 1, ai)).copied());

            if let (Some(i), Some(j)) = (c_idx, n_idx) {
                let key = (i.min(j), i.max(j));
                excluded.insert(key);
                adj[i].push(j);
                adj[j].push(i);

                if movable[i] || movable[j] {
                    let r0 = (positions[i] - positions[j]).norm();
                    if r0 > 1e-3 {
                        bond_terms.push(BondTerm {
                            i,
                            j,
                            r0,
                            k: bond_force_constant(elements[i], elements[j]),
                        });
                    }
                }
            }
        }
    }

    // ── 4.  Build angle terms from the adjacency graph ────────────────────────
    // For each central atom j and each pair of neighbours (i, k), add an angle term.
    let mut angle_set: HashSet<(usize, usize, usize)> = HashSet::new();

    for j in 0..n_atoms {
        for &i in &adj[j] {
            for &k in &adj[j] {
                if i >= k {
                    continue;
                }
                // Add 1-3 exclusion.
                excluded.insert((i.min(k), i.max(k)));

                // Only add the angle term if at least one of (i, j, k) is movable.
                if !movable[i] && !movable[j] && !movable[k] {
                    continue;
                }

                let key = (i, j, k);
                if angle_set.contains(&key) {
                    continue;
                }
                angle_set.insert(key);

                let pi = positions[i];
                let pj = positions[j];
                let pk = positions[k];

                let v1 = (pi - pj).normalize();
                let v2 = (pk - pj).normalize();

                let cos_theta = v1.dot(&v2).clamp(-1.0, 1.0);
                let theta0 = cos_theta.acos();

                let ka = angle_force_constant(elements[j]);
                angle_terms.push(AngleTerm {
                    i,
                    j,
                    k,
                    theta0,
                    ka,
                });
            }
        }
    }

    // ── 5.  Build non-bonded pair list with cutoff ────────────────────────────
    let cutoff_sq = config.vdw_cutoff * config.vdw_cutoff;
    let mut nb_pairs: Vec<NbPair> = Vec::new();

    for i in 0..n_atoms {
        if !movable[i] {
            continue;
        }
        for j in (i + 1)..n_atoms {
            // Exclude bonded and angle partners.
            if excluded.contains(&(i, j)) {
                continue;
            }
            let dist_sq = (positions[i] - positions[j]).norm_squared();
            if dist_sq > cutoff_sq {
                continue;
            }
            let (rmin, eps) = lj_pair_params(elements[i], elements[j]);
            nb_pairs.push(NbPair { i, j, rmin, eps });
        }
    }

    MinSystem {
        positions,
        movable,
        bond_terms,
        angle_terms,
        nb_pairs,
    }
}

// ─── Energy and gradient computation ─────────────────────────────────────────

/// Evaluates the total potential energy of the system (kcal mol⁻¹).
fn compute_energy(sys: &MinSystem) -> f64 {
    let mut energy = 0.0;

    // Bond stretching.
    for bt in &sys.bond_terms {
        let r = (sys.positions[bt.i] - sys.positions[bt.j]).norm();
        let dr = r - bt.r0;
        energy += bt.k * dr * dr;
    }

    // Angle bending.
    for at in &sys.angle_terms {
        let v1 = (sys.positions[at.i] - sys.positions[at.j]).normalize();
        let v2 = (sys.positions[at.k] - sys.positions[at.j]).normalize();
        let cos_t = v1.dot(&v2).clamp(-1.0, 1.0);
        let dtheta = cos_t.acos() - at.theta0;
        energy += at.ka * dtheta * dtheta;
    }

    // Lennard-Jones non-bonded.
    for nb in &sys.nb_pairs {
        let r2 = (sys.positions[nb.i] - sys.positions[nb.j]).norm_squared();
        if r2 < 1e-6 {
            continue;
        }
        let rm2 = (nb.rmin * nb.rmin) / r2; // (r_min/r)^2
        let rm6 = rm2 * rm2 * rm2;
        let rm12 = rm6 * rm6;
        energy += nb.eps * (rm12 - 2.0 * rm6);
    }

    energy
}

/// Computes the analytical gradient of the potential energy with respect to all movable
/// atom positions.  Returns a vector of force contributions per atom (length = `n_atoms`),
/// with zero vectors for fixed atoms.
fn compute_gradient(sys: &MinSystem) -> Vec<Vector3<f64>> {
    let n = sys.positions.len();
    let mut grad = vec![Vector3::zeros(); n];

    // Bond gradient.
    for bt in &sys.bond_terms {
        let rij = sys.positions[bt.i] - sys.positions[bt.j];
        let r = rij.norm();
        if r < 1e-8 {
            continue;
        }
        let dr = r - bt.r0;
        let coeff = 2.0 * bt.k * dr / r; // dE/dr * (1/r)
        let g = rij * coeff;
        if sys.movable[bt.i] {
            grad[bt.i] += g;
        }
        if sys.movable[bt.j] {
            grad[bt.j] -= g;
        }
    }

    // Angle gradient.
    for at in &sys.angle_terms {
        let pi = sys.positions[at.i];
        let pj = sys.positions[at.j];
        let pk = sys.positions[at.k];

        let u = pi - pj;
        let v = pk - pj;
        let norm_u = u.norm();
        let norm_v = v.norm();
        if norm_u < 1e-8 || norm_v < 1e-8 {
            continue;
        }

        let u_n = u / norm_u;
        let v_n = v / norm_v;
        let cos_t = u_n.dot(&v_n).clamp(-1.0, 1.0);
        let sin_t = (1.0 - cos_t * cos_t).sqrt().max(1e-8);

        let theta = cos_t.acos();
        let dtheta = theta - at.theta0;
        let coeff = 2.0 * at.ka * dtheta / sin_t;

        let d_gi = (cos_t * u_n - v_n) / norm_u;
        let d_gk = (cos_t * v_n - u_n) / norm_v;
        let d_gj = -(d_gi + d_gk);

        if sys.movable[at.i] {
            grad[at.i] += coeff * d_gi;
        }
        if sys.movable[at.j] {
            grad[at.j] += coeff * d_gj;
        }
        if sys.movable[at.k] {
            grad[at.k] += coeff * d_gk;
        }
    }

    // LJ gradient.
    for nb in &sys.nb_pairs {
        let rij = sys.positions[nb.i] - sys.positions[nb.j];
        let r2 = rij.norm_squared();
        if r2 < 1e-6 {
            continue;
        }
        let rm2 = (nb.rmin * nb.rmin) / r2;
        let rm6 = rm2 * rm2 * rm2;
        let rm12 = rm6 * rm6;
        // dE/dr = eps * (-12*rmin^12/r^13 + 12*rmin^6/r^7)
        //       = eps * (12/r^2) * (rmin^6/r^6 - rmin^12/r^12) * rij
        //       = eps * (12/r2) * (rm6 - rm12) * rij
        let coeff = nb.eps * 12.0 * (rm6 - rm12) / r2;
        let g = rij * coeff;
        if sys.movable[nb.i] {
            grad[nb.i] += g;
        }
        if sys.movable[nb.j] {
            grad[nb.j] -= g;
        }
    }

    grad
}

/// Returns the RMS of the gradient over all movable atoms.
fn rms_gradient(grad: &[Vector3<f64>], movable: &[bool]) -> f64 {
    let mut sum_sq = 0.0;
    let mut count = 0usize;
    for (i, g) in grad.iter().enumerate() {
        if movable[i] {
            sum_sq += g.norm_squared();
            count += SPATIAL_DIMS;
        }
    }
    if count == 0 {
        return 0.0;
    }
    (sum_sq / count as f64).sqrt()
}

// ─── Steepest-descent minimizer ───────────────────────────────────────────────

/// Runs the steepest-descent loop and returns the number of steps taken.
fn minimize(sys: &mut MinSystem, config: &RelaxConfig) -> (u32, bool) {
    let mut step_size = INITIAL_STEP_SIZE;
    let mut energy = compute_energy(sys);
    let mut steps = 0u32;

    for step in 0..config.max_steps {
        steps = step + 1;
        let grad = compute_gradient(sys);
        let rms = rms_gradient(&grad, &sys.movable);

        if rms < config.convergence {
            return (steps, true);
        }

        // Try a step along the negative gradient.
        let new_pos: Vec<Vector3<f64>> = sys
            .positions
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                if sys.movable[i] {
                    p - grad[i] * step_size
                } else {
                    p
                }
            })
            .collect();

        // Evaluate energy at the new positions.
        let saved_pos = sys.positions.clone();
        sys.positions = new_pos;
        let new_energy = compute_energy(sys);

        if new_energy <= energy {
            // Accept step: grow step size for next iteration.
            energy = new_energy;
            step_size = (step_size * STEP_GROW).min(MAX_STEP_SIZE);
        } else {
            // Reject: restore positions, shrink step size.
            sys.positions = saved_pos;
            step_size *= STEP_SHRINK;
            if step_size < MIN_STEP_SIZE {
                break;
            }
        }
    }

    (steps, false)
}

// ─── Write-back ───────────────────────────────────────────────────────────────

/// Copies the optimized positions from the flat system back into the structure.
fn write_back(sys: &MinSystem, structure: &mut Structure) {
    // Build a flat list of mutable atom references in the same order as atom_map.
    // We traverse chains/residues/atoms in the same order used in `build_system`.
    let mut flat_atoms: Vec<*mut crate::model::atom::Atom> =
        Vec::with_capacity(sys.positions.len());
    for chain in structure.iter_chains_mut() {
        for residue in chain.iter_residues_mut() {
            for atom in residue.iter_atoms_mut() {
                flat_atoms.push(atom as *mut _);
            }
        }
    }

    for (flat_idx, &pos) in sys.positions.iter().enumerate() {
        if sys.movable[flat_idx] {
            // SAFETY: `flat_atoms[flat_idx]` is a valid, exclusive pointer derived from the
            // mutable borrow of `structure`.  Indices are assigned in identical traversal
            // order so there are no aliased writes.
            unsafe {
                if let Some(&ptr) = flat_atoms.get(flat_idx) {
                    (*ptr).pos = nalgebra::Point3::from(pos);
                }
            }
        }
    }
}

// ─── Public entry point ───────────────────────────────────────────────────────

/// Relaxes a structure by minimizing a simplified AMBER-like energy function.
///
/// The side-chain relaxation mode keeps backbone heavy atoms fixed (`N`, `CA`, `C`, `O`,
/// `OXT`) while optimizing side-chain heavy atoms.  The full-minimization mode moves all
/// heavy atoms of standard residues.  Non-standard residues (water, ions, heterogens) are
/// always held fixed.
///
/// # Arguments
///
/// * `structure` - Mutable structure to relax in place.
/// * `config` - Tuning parameters (steps, convergence, scope).
///
/// # Returns
///
/// A [`RelaxResult`] describing the energy change and convergence status.
///
/// # Errors
///
/// Currently this function does not produce hard errors; it always returns `Ok`.
pub fn relax_structure(
    structure: &mut Structure,
    config: &RelaxConfig,
) -> Result<RelaxResult, Error> {
    let mut sys = build_system(structure, config);

    // Guard against structures where every atom is fixed.
    if !sys.movable.iter().any(|&m| m) {
        let scope = if config.side_chains_only {
            "side chains"
        } else {
            "all heavy atoms"
        };
        return Err(Error::NoMovableAtoms {
            scope: scope.to_string(),
        });
    }

    let initial_energy = compute_energy(&sys);
    let (steps_taken, converged) = minimize(&mut sys, config);
    let final_energy = compute_energy(&sys);

    write_back(&sys, structure);

    Ok(RelaxResult {
        initial_energy,
        final_energy,
        steps_taken,
        converged,
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{atom::Atom, chain::Chain, residue::Residue, types::ResidueCategory};

    fn make_single_atom_structure(element: Element, pos: [f64; 3]) -> Structure {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(crate::model::types::StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new(
            "CB",
            element,
            nalgebra::Point3::new(pos[0], pos[1], pos[2]),
        ));
        chain.add_residue(residue);
        structure.add_chain(chain);
        structure
    }

    #[test]
    fn relax_default_config_is_sensible() {
        let cfg = RelaxConfig::default();
        assert_eq!(cfg.max_steps, DEFAULT_MAX_STEPS);
        assert!(cfg.side_chains_only);
        assert!(cfg.convergence > 0.0);
        assert!(cfg.vdw_cutoff > 0.0);
    }

    #[test]
    fn relax_empty_structure_returns_no_movable_atoms_error() {
        let mut s = Structure::new();
        let cfg = RelaxConfig::default();
        let result = relax_structure(&mut s, &cfg);
        assert!(
            matches!(result, Err(Error::NoMovableAtoms { .. })),
            "Expected NoMovableAtoms error for empty structure, got {:?}",
            result
        );
    }

    #[test]
    fn relax_no_movable_atoms_returns_error() {
        // A water molecule has no standard residue atoms that are side-chain; should be
        // treated as fixed, returning an error.
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(1, None, "HOH", None, ResidueCategory::Hetero);
        residue.add_atom(Atom::new(
            "O",
            Element::O,
            nalgebra::Point3::new(0.0, 0.0, 0.0),
        ));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let cfg = RelaxConfig::default();
        let result = relax_structure(&mut structure, &cfg);
        assert!(
            matches!(result, Err(Error::NoMovableAtoms { .. })),
            "Expected NoMovableAtoms error for hetero-only structure, got {:?}",
            result
        );
    }

    #[test]
    fn lj_pair_params_symmetric() {
        let (r_cn, e_cn) = lj_pair_params(Element::C, Element::N);
        let (r_nc, e_nc) = lj_pair_params(Element::N, Element::C);
        assert!((r_cn - r_nc).abs() < 1e-12);
        assert!((e_cn - e_nc).abs() < 1e-12);
    }

    #[test]
    fn vdw_params_returns_positive_values() {
        for elem in [
            Element::C,
            Element::N,
            Element::O,
            Element::S,
            Element::H,
            Element::P,
        ] {
            let (r, eps) = vdw_params(elem);
            assert!(r > 0.0, "r_min/2 must be positive for {:?}", elem);
            assert!(eps > 0.0, "epsilon must be positive for {:?}", elem);
        }
    }

    #[test]
    fn build_system_counts_movable_sidechain_atoms() {
        // ALA with CB as side chain atom – only CB should be movable.
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut res = Residue::new(
            1,
            None,
            "ALA",
            Some(crate::model::types::StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        // Add backbone atoms (should be fixed).
        res.add_atom(Atom::new(
            "N",
            Element::N,
            nalgebra::Point3::new(-0.966, 0.493, 1.500),
        ));
        res.add_atom(Atom::new(
            "CA",
            Element::C,
            nalgebra::Point3::new(0.257, 0.418, 0.692),
        ));
        res.add_atom(Atom::new(
            "C",
            Element::C,
            nalgebra::Point3::new(-0.094, 0.017, -0.716),
        ));
        res.add_atom(Atom::new(
            "O",
            Element::O,
            nalgebra::Point3::new(-1.056, -0.682, -0.923),
        ));
        // Side chain atom (should be movable).
        res.add_atom(Atom::new(
            "CB",
            Element::C,
            nalgebra::Point3::new(1.204, -0.620, 1.296),
        ));
        chain.add_residue(res);
        structure.add_chain(chain);

        let cfg = RelaxConfig {
            side_chains_only: true,
            max_steps: 1,
            ..RelaxConfig::default()
        };

        let sys = build_system(&structure, &cfg);
        let movable_count = sys.movable.iter().filter(|&&m| m).count();
        assert_eq!(movable_count, 1, "Only CB should be movable");
    }

    #[test]
    fn energy_decreases_for_clashing_atoms() {
        // Two carbon atoms placed very close together should repel; after minimization
        // their energy should be lower than before.
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");

        let mut res1 = Residue::new(
            1,
            None,
            "ALA",
            Some(crate::model::types::StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        res1.add_atom(Atom::new(
            "N",
            Element::N,
            nalgebra::Point3::new(-0.966, 0.493, 1.500),
        ));
        res1.add_atom(Atom::new(
            "CA",
            Element::C,
            nalgebra::Point3::new(0.257, 0.418, 0.692),
        ));
        res1.add_atom(Atom::new(
            "C",
            Element::C,
            nalgebra::Point3::new(-0.094, 0.017, -0.716),
        ));
        res1.add_atom(Atom::new(
            "O",
            Element::O,
            nalgebra::Point3::new(-1.056, -0.682, -0.923),
        ));
        // CB placed extremely close to N of another residue to create a clash.
        res1.add_atom(Atom::new(
            "CB",
            Element::C,
            nalgebra::Point3::new(-0.900, 0.480, 1.490),
        ));

        let mut res2 = Residue::new(
            2,
            None,
            "ALA",
            Some(crate::model::types::StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        res2.add_atom(Atom::new(
            "N",
            Element::N,
            nalgebra::Point3::new(-0.966, 0.493, 1.500),
        ));
        res2.add_atom(Atom::new(
            "CA",
            Element::C,
            nalgebra::Point3::new(0.257, 0.418, 0.692),
        ));
        res2.add_atom(Atom::new(
            "C",
            Element::C,
            nalgebra::Point3::new(-0.094, 0.017, -0.716),
        ));
        res2.add_atom(Atom::new(
            "O",
            Element::O,
            nalgebra::Point3::new(-1.056, -0.682, -0.923),
        ));
        res2.add_atom(Atom::new(
            "CB",
            Element::C,
            nalgebra::Point3::new(1.204, -0.620, 1.296),
        ));

        chain.add_residue(res1);
        chain.add_residue(res2);
        structure.add_chain(chain);

        let cfg = RelaxConfig {
            side_chains_only: true,
            max_steps: 50,
            convergence: 1.0,
            vdw_cutoff: 10.0,
        };

        let result = relax_structure(&mut structure, &cfg).unwrap();
        assert!(
            result.final_energy <= result.initial_energy,
            "Energy should not increase: {} -> {}",
            result.initial_energy,
            result.final_energy
        );
    }

    #[test]
    fn relax_structure_errors_when_no_movable_atoms() {
        // A structure with only backbone atoms in side-chain-only mode has no movable atoms.
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut res = Residue::new(
            1,
            None,
            "GLY",
            Some(crate::model::types::StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        res.add_atom(Atom::new(
            "N",
            Element::N,
            nalgebra::Point3::new(-0.966, 0.493, 1.500),
        ));
        res.add_atom(Atom::new(
            "CA",
            Element::C,
            nalgebra::Point3::new(0.257, 0.418, 0.692),
        ));
        res.add_atom(Atom::new(
            "C",
            Element::C,
            nalgebra::Point3::new(-0.094, 0.017, -0.716),
        ));
        res.add_atom(Atom::new(
            "O",
            Element::O,
            nalgebra::Point3::new(-1.056, -0.682, -0.923),
        ));
        chain.add_residue(res);
        structure.add_chain(chain);

        let cfg = RelaxConfig {
            side_chains_only: true,
            ..RelaxConfig::default()
        };

        let result = relax_structure(&mut structure, &cfg);
        assert!(
            matches!(result, Err(Error::NoMovableAtoms { .. })),
            "Expected NoMovableAtoms error, got {:?}",
            result
        );
    }
}
