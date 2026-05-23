/**
 * @file Relax step
 *
 * Pipeline step for coordinate relaxation.
 */

"use client";

import { useShallow } from "zustand/react/shallow";
import { Checkbox, Input } from "@/ui/primitives";
import { ZapIcon } from "@/ui/icons";
import { usePipelineStore } from "@/state";
import { DEFAULT_RELAX_SETTINGS } from "@/core/pipeline";
import { StepWrapper } from "./step-wrapper";

// ============================================================================
// Component
// ============================================================================

export function StepRelax() {
  const { enabled, setEnabled, config, setConfig } = usePipelineStore(
    useShallow((s) => ({
      enabled: s.relaxEnabled,
      setEnabled: s.setRelaxEnabled,
      config: s.relaxConfig,
      setConfig: s.setRelaxConfig,
    }))
  );

  return (
    <StepWrapper
      icon={<ZapIcon className="size-4" />}
      title="Relax"
      enabled={enabled}
      onToggle={setEnabled}
      onReset={() => setConfig(DEFAULT_RELAX_SETTINGS)}
    >
      <Checkbox
        label="Side chains only"
        description="When disabled, all standard-residue heavy atoms are relaxed."
        checked={config.sideChainsOnly}
        onChange={(e) => setConfig({ sideChainsOnly: e.target.checked })}
      />

      <div className="grid grid-cols-2 gap-3">
        <Input
          label="Max steps"
          type="number"
          min={1}
          step={1}
          value={config.maxSteps}
          onChange={(e) =>
            setConfig({ maxSteps: Math.max(1, parseInt(e.target.value) || 200) })
          }
        />
        <Input
          label="Convergence"
          type="number"
          min={0}
          step={0.1}
          value={config.convergence}
          onChange={(e) =>
            setConfig({ convergence: Math.max(0, parseFloat(e.target.value) || 1) })
          }
        />
      </div>

      <Input
        label="VDW cutoff (Å)"
        type="number"
        min={0}
        step={0.1}
        value={config.vdwCutoff}
        onChange={(e) =>
          setConfig({ vdwCutoff: Math.max(0, parseFloat(e.target.value) || 10) })
        }
      />
    </StepWrapper>
  );
}
