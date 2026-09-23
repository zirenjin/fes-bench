# Hf hcp QH gate audit

The finite-displacement diagnosis is in
`results/phase2/hcp_imaginary_diagnosis.md`.  The SSE-PBE head has scattered
negative modes (101/133 negative q-points; 0 of 162 negative modes in the
5%-radius Γ neighborhood), so ASR and Γ exclusion were not used.  The
metallic-domain `Domains_Alloy` head has no mode below -0.05 THz (minimum
0.000777913 THz) and is the selected hcp representative head.

For the fresh Domains_Alloy QH rerun, the hcp `F_QH-reference` residual was
evaluated on all 1,252 shared temperatures.  Its range is
`[-7.1854488667, -7.0853388116] eV/atom`; the large nearly constant offset is kept
as a documented model/reference offset rather than hidden.  The temperature
variation is smooth: maximum absolute first difference
`2.0824305e-4 eV/atom` per grid step and maximum absolute second difference
`1.0023968e-4 eV/atom^2`.
There are no positive-step excursions in the bcc comparison (`0`), while hcp
has 73 tiny downward steps over the tabulated grid; these do not create a
non-smooth jump or branch discontinuity.

**Gate decision:** `hcp qh_reliable=true` for the selected Domains_Alloy head,
subject to retaining the diagnosis, single-point sanity output, and rerun
artifacts.  `bcc qh_reliable=false` remains unchanged.  This gate does not
claim absolute agreement with the reference; it only records the requested
absence of scattered imaginary modes and a smooth temperature deviation.
