# Next — abjm-c123

*Updated 2026-10-11 (evening)*

**Focus:** tree-level structure constants in ABJM theory, full su(4) sector. Active thread: the su(4) alternating
spin chain (fund ⊗ anti-fund ⊗ …) — SoV basis, localisation, free-parameter rotations on few sites — now being
connected to the ABJM 3-pt function of the note `StructuringABJM.pdf` (repo root) following arXiv:2603.26617.
Everything is in `Paul/Mathematica/su4-alternating/`: **read its `README.md`, then `ABJM-3pt-transfer-matrices.md`**.

## Where it stands

- L = 2: SoV basis ⟨0|t_a(θ₁−h)t₁(θ₂)^k, gl(4) B-operator, qKZ localisation — done, exact (notebook).
- L = 4, companion frame: one free κ per site on the bare SoV vacuum; {1,3} common rotation = global exp(E₁₂);
  good-reference criterion w G vᵀ = 0; two commuting common steps; 3⊗3 word span; the 2-parameter state only allows
  end-growing (fixed γ) rotations. (README, numpy.)
- **ABJM 3-pt (diagonal frame, L = 2 pairs):** the note's (3.6) = ⟨W₋₁^L| with the first pairs rotated by 1 − Λ/z₄,
  exactly 2603's O₂ = O₃·(1 − κ·twist) on a block. Odd bridge (one field) works exactly:
  C ∝ (1 − z₄/(c₁τ₁(θ₁)))/(1 − z₄/z₁)·⟨W₋₁^L|Ψ⟩, verified on on-shell HW states through the untwisting limit.
  Explicit L = 2 table (all bridges): adjoints [1,0,1] only seen by odd bridges (C₁ = 1/√8, 1/√24), singlets only by even. Even bridge (a pair) is **not** reached by localising words — missing piece: Λ on the
  F site alone. Paul's HW trick (free global e^{bE₁₄}) gives an exact finite-twist one-word formula with
  b = z₁/(z₁−z₄), but b diverges and the untwisting limit is wrong (state-dependent F/T); 2603 avoids this with a
  continuous κ-family, which for ABJM pairs is exactly what is missing.

## Next steps

1. Find the pair κ-family (F-site-only twist rotation at the left end, or a combination cancelling the singlet of
   ⟨W|t₁(θ₂)); try L = 3 pairs. Then K_κ for ℓ pairs and untwisting with κ → 1.
2. **Resume here:** odd bridges at L = 3 pairs — `abjm-bridges.wb` §3 with `Lp = 3`, **at the ABJM point
   θ_A − θ_F = 2h** (found 2026-10-12: only there H_MZ commutes with the t's; θ = 0, used in all earlier bridge
   computations, is singular — S drops out of the t-algebra; L = 2 table redone, see ABJM-md §5.7 box). Also check
   whether earlier ABJM results (§5.5–5.6, untwist scripts) assumed θ = 0. Done 2026-10-12: t₃ in vector form
   (`dual-lax-t3.wb`: ℒ̄_A = (u−θ+h) − hP, ℒ̄_F = (u−θ−3h) + hK, twist det G·G^{−T}; 0.16 s per pass at L = 6 sites), t₂ in
   vector form (same notebook §4, ∧² Lax, F shifted by h; `t2V` in `abjm-bridges.wb` — add it to `mixT` at L = 3), and
   the L = 2 degenerate matching (needs momentum: pair-periodic θ + twisted shift Ŝ = S∘ρ_{3..L}(Λ); t₁, t₂, t₃ do not
   separate the 1/√8 state from an odd-momentum partner). See `ABJM-3pt-transfer-matrices.md` §5.7.
2b. **LR symmetry (2026-10-12, `abjm-LR-symmetry.wb`):** all HW cyclic states up to L = 4 pairs saved in
   `su4-alternating/data/abjm-hw-states.wxf` — load, don't recompute. No 2603-type LR selection rule for ⟨B_ℓ|;
   vanishing = R-charge ([p,0,p], |ℓ−L| ≤ p) + bridge parity (L+p odd ⇒ ℓ ≡ L+p mod 2) + two isolated zeros at L = 4
   (E = 32 [1,0,1], E = 8 [2,0,2]). Open: 𝒫-symmetrised bra (both bridge placements)? twisted version?
3. Integrability of ⟨W_a^L| (note's (4.8)); later the ⟨W^L|Ψ⟩ determinant (Q stage).
4. Notebook versions of the numpy results ({1,3} common rotation closed form, frame map, odd-bridge formula).

## Housekeeping

- The two `.wb` notebooks have no native outputs — rerun in VS Code to store them.
- Committed up to 56b9407 (thread + `chk_rotref.py`); the ABJM write-up, README/NEXT updates and the new `checks/`
  scripts (`rotbase.py`, `chk_span13_words.py`, `scan_*.py`, `chk_frame.py`, `abjm_*.py`) are **not committed yet**
  (all smoke-tested, exit 0).
- 2603 LaTeX source is not in the repo (read from arXiv e-print 2603.26617v1, `short.tex`).
