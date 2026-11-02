# Thermodynamic map — comparison and aggregation policy

This map combines navigation over different experiments, not their statistical scores. Each contract retains its reference, task, sample and input budget. The active overview is one task-by-method matrix. Percentage accuracy uses a 0–100 colour scale. Error colours normalize only within each task row; diagnostic values use neutral colours. MAE, RMS, samples and reference levels remain distinct. Older detailed contract records remain downloadable but no longer form separate displayed heatmaps. No overall score is computed.

## Aggregation correction

Checking `analyze_covariance.py` against the exported JSON shows that Alexandria Q1 single-phase RMS uses an unweighted average over entries; Q1 pair RMS and Q2–Q4 RMS use an unweighted average over pairs. They are phase-weighted and pair-weighted respectively. Cancellation and Q3 ordering summaries use `boot()` with formula-level means. The former blanket description of all audit metrics as formula-weighted was too broad. The November 1 page captions and source table now distinguish these quantities. No scientific values were changed.

[Audit aggregation implementation](../updates-20261101/sources/aggregation-analysis.py) · [Original result fields](../updates-20261101/sources/alexandria-primary.json)

## Coverage and missingness

Workflow coverage uses all intended pairs (295/351 for the relaxed harmonic pipeline). Selective-decision coverage uses eligible non-degenerate pairs under the archived abstention rule. They are not interchangeable. F1 imaginary-mode counts do not certify its end-to-end workflow coverage; that field remains NR. F1b remains PENDING as of this frozen website snapshot; no remote training status was queried.

NR = not reported or not run; RF = independent reference missing; NA = outside the task; NE = no eligible matched event; PENDING = awaiting result. The legacy map retains its own older state definitions, including NE for no reference event.

## Physics and inference

Carbon labels are single-phase harmonic Svib at fixed volume and 300 K. SiO2 slope errors are effective relative temperature response from G curves. Neither is relabelled as the other. Alexandria fixed-volume Helmholtz F is not Gibbs G. SiO2 source-arm comparisons are development results with unequal source supervision, not an untouched model selection test. Frequency counterfactuals are ideal spectral error removals and not predicted gains from training.

All values come from the linked versioned source records. The machine-readable JSON and CSV provide per-cell provenance; no new DFT, MD, phonon inference or model training is performed by the website.
