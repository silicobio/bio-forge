use anyhow::{Context, Result};
use clap::Args;

use bio_forge::Structure;
use bio_forge::ops::{RelaxConfig, relax_structure};

use crate::commands::run_with_spinner;

/// Minimizes side-chain or whole-structure energy using a simplified AMBER-like force field.
///
/// By default only protein side-chain heavy atoms are moved while backbone atoms (`N`, `CA`,
/// `C`, `O`, `OXT`) and all non-standard residues (water, ions, ligands) remain fixed.
/// Pass `--full` to also optimize backbone heavy atoms of standard residues.
#[derive(Debug, Args)]
pub struct RelaxArgs {
    /// Move all standard-residue heavy atoms instead of side chains only.
    #[arg(long = "full", default_value_t = false)]
    pub full: bool,

    /// Maximum number of steepest-descent minimization steps.
    #[arg(long = "steps", default_value_t = 200)]
    pub steps: u32,

    /// RMS-gradient convergence threshold (kcal mol⁻¹ Å⁻¹).
    #[arg(long = "convergence", default_value_t = 1.0)]
    pub convergence: f64,

    /// Lennard-Jones non-bonded cutoff distance (Å).
    #[arg(long = "vdw-cutoff", default_value_t = 10.0)]
    pub vdw_cutoff: f64,
}

/// Runs the relaxation pipeline and prints a short summary to stderr.
pub fn run(structure: &mut Structure, args: &RelaxArgs) -> Result<()> {
    let config = RelaxConfig {
        max_steps: args.steps,
        side_chains_only: !args.full,
        convergence: args.convergence,
        vdw_cutoff: args.vdw_cutoff,
    };

    let scope = if config.side_chains_only { "side chains" } else { "full structure" };
    let message = format!("Relaxing {scope} ({} steps max)", config.max_steps);

    let result = run_with_spinner(&message, || {
        relax_structure(structure, &config).context("Failed to relax structure")
    })?;

    eprintln!(
        "  Energy: {:.2} → {:.2} kcal/mol  ({} steps{})",
        result.initial_energy,
        result.final_energy,
        result.steps_taken,
        if result.converged { ", converged" } else { "" }
    );

    Ok(())
}
