/**
 * @file WASM type definitions
 *
 * TypeScript interfaces for the BioForge WASM bindings.
 */

export type {
  // Data structures
  StructureInfo,
  ChainInfo,
  // Configuration types
  CleanConfig,
  HydroConfig,
  RelaxConfig,
  SolvateConfig,
  TopologyConfig,
  // Classes
  Structure as WasmStructure,
  Template as WasmTemplate,
  Topology as WasmTopology,
} from "bio-forge";

/** Chain polymer classification */
export type PolymerType = "protein" | "nucleic" | "solvent" | "hetero";

/** Histidine tautomer selection strategy */
export type HisStrategy = "hid" | "hie" | "random" | "network";

/** Cation species for solvation */
export type CationSpecies = "Na" | "K" | "Mg" | "Ca" | "Li" | "Zn";

/** Anion species for solvation */
export type AnionSpecies = "Cl" | "Br" | "I" | "F";

import type { Structure, Template, Topology } from "bio-forge";

/** The complete WASM module interface */
export interface WasmModule {
  Structure: typeof Structure;
  Topology: typeof Topology;
  Template: typeof Template;
}
