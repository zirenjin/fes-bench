# FES residual-head descriptor provenance

## Conclusion

The residual correction input `z` is the descriptor tensor produced by the
shared DPA descriptor backbone, not a head-specific output.  The FES baseline
and residual correction consume the same tensor in one atomic-model forward.
The selected DPA multitask head controls the fitting branch (`E` baseline and
the FES correction attached to that branch); it does not redefine `z`.

## Code-path evidence

The evidence below is from the `feat/fes-head` checkout at
`/home/ziren/aisi-intern/deepmd-kit`:

1. `deepmd/pt/utils/multi_task.py`, in the shared-item handling of
   `MultiTaskEmbedding`, treats `descriptor` as a supported
   shared item, expands each branch from `shared_dict`, and records links with
   `shared_type: "descriptor"` and a `shared_level`.  Thus a multitask model
   can share descriptor parameters while retaining separate fitting branches.
2. `deepmd/pt/model/atomic_model/dp_atomic_model.py`, in the atomic-model
   `forward` path, calls
   `self.descriptor(...)` once and stores its result in the local tensor
   `descriptor`; the same method then calls `self.fitting_net(...)` with that
   tensor.
3. `deepmd/pt/model/task/free_energy.py`, in `FreeEnergyFittingNet.forward`,
   passes `descriptor` to the frozen `baseline` and initializes
   `corr_descriptor = descriptor`
   (optionally concatenating only the encoded thermodynamic state) before
   passing it to `correction`.  There is no second head-specific descriptor
   extraction.
4. `deepmd/pt/model/model/make_model.py`, in the model-forward assembly,
   invokes
   `forward_common_atomic` once per model forward, so the two FES subnets are
   downstream consumers of the same descriptor pass.
5. `deepmd/pt/model/model/free_energy_model.py`, in its state-vector
   construction and forward methods, constructs the full `[T, P, v, c]`
   state vector separately and forwards it through
   the ordinary FES model path.  `F_QH` therefore enters the correction state,
   not the descriptor definition.

## Fixed method configuration

`configs/models/head_policy.yaml` fixes the method statement as:

> `z` is the shared DPA-3.1-3M descriptor trunk output; `E` and `F_QH` are
> evaluated through the system-selected head.

The recorded choices are `Domains_Alloy` for both `energy_head` and `fqh_head`
for Hf, SiO2, Ti, and Zr.  Hf is supported by the completed hcp reciprocal-
space diagnosis; Ti/Zr are metallic-domain choices whose representative
sanity/QH results are recorded in their Phase 1/2 audits.  The Hf bcc hold is
retained, and no architecture change was made.

## Deviations from design

- The checkpoint's serialized branch names and the FES code's `shared_dict`
  links are documented separately: selecting a head selects a fitting branch;
  it is not evidence for a private descriptor latent.
- No new descriptor or residual architecture was introduced.  This deliverable
  records the existing code path and freezes the selection policy only.
