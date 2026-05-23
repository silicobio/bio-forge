/**
 * @file WASM module exports
 */

export { initWasm, isWasmReady, getWasm, requireWasm } from "./module";

export {
  getStructureInfo,
  cleanStructure,
  repairStructure,
  addHydrogens,
  relaxStructure,
  solvateStructure,
  buildTopology,
} from "./structure";
export type {
  StructureFormat,
  WasmStructure,
  WasmTopology,
  WasmTemplate,
} from "./structure";

export { exportTopology, getBondCount } from "./topology";
export type { TopologyConfig } from "./topology";

export type {
  WasmModule,
  StructureInfo,
  ChainInfo,
  PolymerType,
  CleanConfig,
  HydroConfig,
  HisStrategy,
  RelaxConfig,
  SolvateConfig,
  CationSpecies,
  AnionSpecies,
} from "./types";
