# Narrative revision: findings first, diagnostics second

## Canonical terms

| Term | Meaning |
| --- | --- |
| Source fit | Error on permitted source training states |
| Transfer loss | Increase in error in the declared held-out task |
| Information loss | Thermodynamic quantity not preserved by this prediction |
| Optimization Sanity Check | Source-fit and held-out error audit; not a convergence certificate |

## Claim–evidence map

1. **Source fit and transfer are distinct.** C-DPA fits source gaps closely, while held-out target-pair MAE rises in both reported folds. Supported by frozen C point predictions and saved readout matrices.
2. **Mean ΔG MAE and phase decisions can rank methods differently.** The SiO₂ Phase-ID/T3 comparison retains both seeds. Supported by existing paired score rows.

## Result allocation

- Core findings: opening two cards, C source/target table and the paired phase-decision example.
- Necessary support: source-fit empirical distributions and training/validation/test summary statistics.
- Material qualification: the failed legacy harmonic baseline remains explicit beside the Hf result.
- Scope/statistical units: consolidated in one expandable section.
- Mechanism and layer choice: proposed probe protocol; no new result is claimed.
- Provenance and historical audits: remain in Methods, figure source data and supporting reports.

## Replaced or compressed

The opening causal disclaimer is replaced by two measured findings. The phrases “real counterexample, with a narrow scope” and “one fold ... does not establish broad material performance” are removed from the generated main narrative. Quantitative task definitions and the invalid-baseline warning remain. The new source-fit module replaces repeated cautions with observable distributions.

## Text budget

| Section | Before words | After words |
| --- | ---: | ---: |
| overview | 117 | 223 |
| failure-landscape | 676 | 642 |
| information-loss | 546 | 513 |
| baseline-contribution | 547 | 701 |
| next | 189 | 354 |

The source-fit section grows because it adds previously missing evidence. The overview now states the findings before their scope.

## Experiment status

Existing predictions were reanalysed; no model was retrained. Intermediate-layer extraction and probes remain unrun. Saved C source contrast matrices have ranks 1 and 2 for DPA; the metal C-simple contrast rank is 1 despite two input contrast vectors. Readout hashes are checked against the saved model-selection record.

## Update after the P0 transfer audit

The opening now separates three result states: observed known-phase fits, observed MAE/phase-decision disagreement, and audited but preliminary transfer. The two main findings are metric disagreement and the coexistence of high envelope accuracy with wrong pair-crossing topology. The new target phase changes 917/1649 decisions; the source-only explanation is rejected. A real Ti family-control phase-order bug is repaired in the derived V2 score layer. The one-contrast ridge rank constraint now precedes the layer-probing plan. Source-fit distributions remain necessary support.
