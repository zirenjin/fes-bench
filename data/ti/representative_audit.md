# Ti representative audit

Both phases use Materials Project structures from the same elemental Ti source
family: hcp `mp-46` and bcc `mp-73`.  The MP records identify hcp as
`P6_3/mmc` (No. 194) with `a=2.935 Å, c=4.659 Å`, and bcc as `Im-3m`
(No. 229) with `a=3.245 Å`.  The source extxyz files are deterministic
ideal-prototype realizations of those MP entries (using the recorded MP lattice
constants and fractional prototype positions), with no additional atoms or
distortions.

| phase | source | space group before relaxation | space group after relaxation | head/protocol |
|---|---|---|---|---|
| hcp | Materials Project `mp-46` | P6₃/mmc (194) | P6₃/mmc (194) | Domains_Alloy; FIRE, cell relaxed, converged |
| bcc | Materials Project `mp-73` | Im-3m (229) | Im-3m (229) | Domains_Alloy; FIRE, cell relaxed, converged |

The hcp relaxation changed the energy from `-6.43146053` to `-6.43206096`
eV/atom and volume from `17.37839073` to `17.55476253 Å³/atom`; the bcc
values changed from `-6.33666357` to `-6.33831224` eV/atom and
`17.08496556` to `17.38597702 Å³/atom`.

The thermodynamic reference remains DaRUS `10.18419/DARUS-3582`; it is not
being replaced by the MP structure source.

The isolated Domains_Alloy QH diagnostic reports `qh_reliable=true` for hcp
and `qh_reliable=false` for bcc because the latter has the expected imaginary
branches.  Details are in `hcp/phonon_report.json` and `bcc/phonon_report.json`.
