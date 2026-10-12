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
2c. **Q-functions + configurations (2026-10-12):** `abjm-Qfunctions.wb` — Q₁, Q₁₂, Q₁₂₃ for all 184 states
   (`data/abjm-Qfunctions.wxf`); LR ⇔ Q₁(u) ∝ Q₁₂₃(u+h), spectrum = LR singlets + LR doublets (as in 2603).
   `abjm-3pt-configs.wb` — no two-block BPS configuration (any polarisations, ℓ, A/F start) gives an LR selection
   rule at L = 3, 4; "only non-LR couple" happens in many configurations but not uniformly in L. Next idea: bras
   annihilated by t̄₁(u) − t̄₃(u+h) (twisted, à la 2603 ω = 1), or bridges built from transfer matrices.
2d. **Factorised LR-annihilated bras (`abjm-LR-factorized.wb`, L = 2):** ⟨B|(t̄₁(u) − t̄₃(u+h)) = 0 has a large
   kernel (136/256 untwisted; 36/256 at twist z₁z₄ = z₂z₃ = 1), but the only product states in it are trivial
   eigen-bras (untwisted: BPS orbit; twisted: e_a ⊗ e_{5−a} vacua); the note's W_a only at a = 0. Next: L = 3,
   non-diagonal LR-compatible twists.
2e. **Other discrete symmetries (`abjm-other-symmetries.wb`):** besides LR (𝒫): parity P (𝒬): τ̄₁(u) = τ̄₁(2h−u)
   ⇔ Q₁(u) ∝ Q₁₂₃(2h−u); LP (𝒫𝒬): τ̄₁(u) = τ̄₃(3h−u) ⇔ Q₁(u) ∝ Q₁(h−u). All occur (L=4: 24 LR, 60 P, 40 LP of
   152). Product states: P is preserved by EVERY two-site state (v⊗w)^{⊗L} (L = 2, 3, 4, exact to 1e-16) —
   integrable-boundary-state-like; LR none (L=2); LP only an L=2 accident (SU(2)xSU(2) sector). (Per-length
   product searches abandoned: the 3-pt bras must use the SAME polarisations at every L — see 2f.)
2f. **Length-independent 3-pt configurations with a selection rule (`abjm-3pt-configs.wb`, last cells):**
   P-selection rules (C = 0 unless Q₁(u) ∝ Q₁₂₃(2h−u)), operator-level on the cyclic sector at L = 3, 4:
   note's configuration with bridge ℓ = 1 (or 2L−1); half-degenerate orbits at ℓ = 1; both-pairings-zero orbit
   (block 2 = (e₂,e₃)) at ℓ = 2 (or 2L−2). No uniform LR or LP rule in the two-block BPS class.
2g. **SoV at 4 sites (`sov-companion-L4.wb`, notes §8):** companion-twist word basis
   ⟨0|t_{a1}(θ1−h)t1(θ2)^{k2}t_{a3}(θ3−h)t1(θ4)^{k4} has rank 256. Sparse BPS expansions: F-only rotations
   (v = e1, w ∈ e1^⊥) ≤ 9 terms (exp κE34: 4); A-only exp κE12: 13; generic 256. Next: closed forms, map to the
   physical 3-pt polarisations, two-block bras from F-only blocks. B-check: at generic θ all 256 words are
   B-eigen (genuine SoV). At the ABJM point (θ_A,n = θ_F,n + 2h) the RV B is DEGENERATE (169 = 13² distinct
   eigenvalues for d = θ_A − θ_F ∈ {h, 2h, 3h}; simple for d ∈ {0, 4h, generic}); shifted word bases exist but are
   not SoV. Next: limit d → 2h, or a different B (4̄ aux / other frame) that stays simple at d = 2h.
   Write-up of all results: `notes/abjm-transfer-matrices-symmetries.tex/.pdf`.
3. Integrability of ⟨W_a^L| (note's (4.8)); later the ⟨W^L|Ψ⟩ determinant (Q stage).
4. Notebook versions of the numpy results ({1,3} common rotation closed form, frame map, odd-bridge formula).

## Housekeeping

- The two `.wb` notebooks have no native outputs — rerun in VS Code to store them.
- Committed up to 56b9407 (thread + `chk_rotref.py`); the ABJM write-up, README/NEXT updates and the new `checks/`
  scripts (`rotbase.py`, `chk_span13_words.py`, `scan_*.py`, `chk_frame.py`, `abjm_*.py`) are **not committed yet**
  (all smoke-tested, exit 0).
- 2603 LaTeX source is not in the repo (read from arXiv e-print 2603.26617v1, `short.tex`).
