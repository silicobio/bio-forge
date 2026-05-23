/**
 * @file Pipeline panel
 *
 * Main pipeline configuration panel containing all processing steps.
 */

"use client";

import { SparklesIcon } from "@/ui/icons";
import { usePipelineStore } from "@/state";
import { StepClean } from "./step-clean";
import { StepRepair } from "./step-repair";
import { StepRelax } from "./step-relax";
import { StepHydro } from "./step-hydro";
import { StepSolvate } from "./step-solvate";
import { StepTopology } from "./step-topology";

// ============================================================================
// Component
// ============================================================================

export function PipelinePanel() {
  const hasAnyEnabled = usePipelineStore(
    (s) =>
      s.cleanEnabled ||
      s.repairEnabled ||
      s.relaxEnabled ||
      s.hydroEnabled ||
      s.solvateEnabled ||
      s.topologyEnabled
  );

  return (
    <div className="h-full flex flex-col overflow-hidden">
      {/* Header */}
      <div className="p-4 border-b border-border shrink-0">
        <h2 className="text-lg font-semibold flex items-center gap-2">
          <SparklesIcon className="size-5 text-primary" />
          Pipeline
        </h2>
        <p className="text-sm text-muted-foreground mt-1">
          Configure processing steps
        </p>
        {!hasAnyEnabled && (
          <p className="text-xs text-warning mt-2">
            Enable at least one step to process files
          </p>
        )}
      </div>

      {/* Steps */}
      <div className="flex-1 overflow-y-auto p-4 space-y-3">
        <StepClean />
        <StepRepair />
        <StepRelax />
        <StepHydro />
        <StepSolvate />
        <StepTopology />
      </div>
    </div>
  );
}
