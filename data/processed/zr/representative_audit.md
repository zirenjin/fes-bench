# Zr representative audit

Both phases use Materials Project structures from the same elemental Zr source
family: hcp `mp-131` and bcc `mp-41`.  The MP records identify hcp as
`P6_3/mmc` (No. 194) with `a=3.239 Å, c=5.172 Å`, and bcc as `Im-3m`
(No. 229) with `a=3.565 Å`.  The source extxyz files are deterministic
ideal-prototype realizations of those MP entries (using the recorded MP lattice
constants and fractional prototype positions), with no additional atoms or
distortions.

| phase | source | space group before relaxation | space group after relaxation | head/protocol |
|---|---|---|---|---|
| hcp | Materials Project `mp-131` | P6₃/mmc (194) | P6₃/mmc (194) | Domains_Alloy; FIRE, cell relaxed, converged |
| bcc | Materials Project `mp-41` | Im-3m (229) | Im-3m (229) | Domains_Alloy; FIRE, cell relaxed, converged |

The hcp relaxation changed the energy from `-7.15151864` to `-7.15249245`
eV/atom and volume from `23.49530290` to `23.36903808 Å³/atom`; the bcc
values changed from `-7.08675297` to `-7.08700560` eV/atom and
`22.65419356` to `22.79645746 Å³/atom`.

The thermodynamic reference remains DaRUS `10.18419/DARUS-3582`; it is not
being replaced by the MP structure source.

The isolated Domains_Alloy QH diagnostic reports `qh_reliable=true` for hcp
and `qh_reliable=false` for bcc because the latter has the expected imaginary
branches.  Details are in `hcp/phonon_report.json` and `bcc/phonon_report.json`.
