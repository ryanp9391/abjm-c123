# Next — abjm-c123

*Updated 2026-10-11*

**Focus:** tree-level structure constants in ABJM theory, full su(4) sector. Active thread: the su(4) alternating
spin chain (fund ⊗ anti-fund ⊗ …) — SoV basis, localisation, free-parameter rotations on few sites.
Everything is in `Paul/Mathematica/su4-alternating/` — **read its `README.md`** (moved from gaudin-cg 2026-10-11;
nothing there is committed yet).

## Where it stands

- L = 2: SoV basis ⟨0|t_a(θ₁−h)t₁(θ₂)^k (plain transfer matrices), gl(4) B-operator, qKZ localisation — done, exact.
- L = 4 on the bare SoV vacuum (companion twist): every site has one localising transfer matrix → one free κ per site
  (site vectors e₁ᵀ(1−κ₁G⁻ᵀ), e₄ᵀ(1−κ₂G), e₁ᵀ(1−κ₃Gᵀ), e₄ᵀ(1−κ₄G⁻¹)), overlaps are explicit polynomials in
  transfer-matrix eigenvalues × ⟨0|Ψ⟩; global torus ρ(X), X ∈ C[G_cᵀ], adds 3 more (7 parameters total).
- Newest (numpy only): ⟨0|(1 + a t₁(θ₁) + b t₁(θ₃) + c t₁(θ₁)t₁(θ₃)) = ⟨0|ρ₁(g)ρ₃(g), g = 1 − κG_cᵀ — a common rotation
  of the two anti-fundamental sites with a free parameter (`checks/chk_common13.py`).

- Newer (numpy, `checks/chk_rotref.py`): the {1,3} rotation is a global rotation exp(E₁₂) (trivial on F vacua). Rotated
  references v⊗w⊗v⊗w work like ⟨0| iff w G vᵀ = 0; from ⟨0| this allows two commuting common steps, e.g. ⟨0|P₁₃P₂₄ = one global
  rotation on all four sites with factorised overlap. See the README section "Globally rotated reference states".

## Next steps

1. Put the {1,3} common rotation in a notebook exactly (closed form of the constant m).
2. Same for {2,4} (fundamental sites), then try one g on all four sites (resolvent on the F sites).
3. More than one parameter per site; L = 4 SoV basis; Q-function form of the overlaps.
4. Connect to the ABJM structure constants.

## Housekeeping

- The two `.wb` notebooks have no native outputs — rerun in VS Code to store them.
- Nothing from the su4-alternating thread is committed yet.
- The `SU3_Tests` notebooks were already in the repo before this thread (committed in "Sync work in progress").
