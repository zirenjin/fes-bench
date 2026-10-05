# T2 static-term audit

Status: **invalid_existing_results**.

The physical channel is canonical `F_QH_eV_per_atom`, and every phase audit confirms it contains the E_DPA static term. The existing T2 plan declares an additive `E_DPA + F_QH` baseline, so those runs are invalid.

Required correction: `F_QH + r_theta(z, T) + c_system`; retrain tlog seeds 11/23/37.
