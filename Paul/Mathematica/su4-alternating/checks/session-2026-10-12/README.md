# Python backups of the 2026-10-12 agent computations

The Wolfram Engine licence expired mid-session and the shared kernel died, so the `.wb` notebooks written that day
(`twisted-vacuum-*.wb`, `abjm-nongood-reference-limit.wb`, `pair-rotation-*.wb`, `sov-*.wb`, `odd-bridge-general-L.wb`,
`twoblock-polarisation-forms.wb`, `pair-group-rotation.wb`, `even-bridge-*.wb`) have unrun cells with python results pasted as
markdown. These folders are verbatim copies of the agents' python (numpy / fractions / sympy / mpmath) that produced those
results. Rerun the notebooks in VS Code once the licence is renewed; these scripts are the fallback / cross-check.

| folder | notebook | content |
|---|---|---|
| `pairgroup/` | `pair-group-rotation.wb` | `pgbase.py` (conventions, fused t_k, K_κ word), RESULTS.md: the exact 5-word common pair rotation at the ABJM point, z₂ = z₃ |
| `sovq/` | `sov-Q-dictionary-L4.wb` | `sov_l4.py`: 4-site SoV basis, twisted Q's, dictionary, dual ket basis, measure |
| `twoblock/` | `sov-twoblock-expansion-L4.wb` | two-block bras in the SoV basis, ⟨0|_comp ∈ HW orbit of vac^L, bordered Casoratian |
| `onepair/` | `sov-onepair-determinant.wb` | exact one-pair determinant formula, two-pair non-factorisation, F–F exchange formula |
| `oddbridge/` | `odd-bridge-general-L.wb` | general-L chain constants, one-field bridges at L = 3, untwisting, `wxf.py` parser |
| (see `checks/twoblock-polforms/`) | `twoblock-polarisation-forms.wb` | Borel-invariant forms of the two-block overlaps for all saved HW states |
| `locality/` | `pair-rotation-locality-L3.wb` | brute-force locality test at 3 pairs (in progress at copy time) |
| (see `checks/abjmlimit/`) | `sov-abjm-point-limit.wb` | ABJM-point limit of the SoV basis: regularised (u-derivative) basis, exact identities t₃(u₀) = 0, t₂(u₀) = t₁(u₀)t₁(θ_F) |

Scripts write only to the session temp directory; paths inside them may need adjusting.
