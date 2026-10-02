# Invalid T3 polynomial/SiO2/seed-11 result

These raw files were generated before the checkpoint type-map and FES label
scale fixes.  The run used compact local `type.raw` IDs with a full
DPA-3.1-3M periodic-table branch and fed eV/atom labels to a loss that divides
by atom count.  Its 79.35 meV/atom training ΔG MAE is invalid for the current
protocol.  The files are retained for audit only and are not read by table
builders.
