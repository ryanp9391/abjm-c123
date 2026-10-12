# su(4) alternating spin chain — SoV basis, localisation, rotations on few sites

Thread started 2026-10-11 in the gaudin-cg repo and moved here the same day; it is input for the ABJM su(4)-sector
structure constants. Read this first; the notebooks have the details. Evidence tags: **[exact]** (exact rational/symbolic
in a notebook), **[numeric]** (60-digit, or numpy 10⁻¹⁴ where marked), **[numpy only]** (checked in `checks/`, not yet in a
notebook).

## Question (Paul)

Alternating chain fundamental ⊗ anti-fundamental ⊗ fund ⊗ anti-fund of su(4), starting at L = 2. Build the SoV basis in the
companion twist frame (Ryan–Volin), study SoV localisation, and construct rotations that act on only a few sites and
depend on free parameters (as in Bargheer–Bercini–Lefundes–Ryan, arXiv:2603.26617, App. F–G, M). **Finite generic twist
only — no untwisting/ε in this thread** (Paul).

## Conventions and naming

- Lax ℒ⁽ˢ⁾(u) = Σ e_ij ⊗ ((u−θ_s)δ_ij + h ρ_s(E_ij)) (gaudin-cg "new conventions"), T = ℒ⁽¹⁾ℒ⁽²⁾⋯ (site 1 leftmost in the
  auxiliary product), twist on the right T^G = T·G, t_a(u) = Σ_{|I|=a} qmin_{I,I}[T^G](u) (rows antisymmetrised, column k
  at u−(k−1)h). Companion twist G_c(z) = {{c₁,−c₂,c₃,−c₄},{1,0,0,0},{0,1,0,0},{0,0,1,0}}, c_k = e_k(z).
- **Paul's naming:** a *fundamental* site is the one whose SoV layer is built with t₁ⁿ (n = 0..3): ρ_F(E_ij) = −e_ji, HW e₄,
  group acts by g^{−T}. An *anti-fundamental* site is built with t_n (n = 0..3): ρ_A(E_ij) = e_ij, HW e₁, group acts by g.
  (This is the λ ↔ λ* swap of the new Lax: u + hΣe_ij⊗ρ(E_ij) is Ryan–Volin's u − ħΣE_ij⊗π(E_ji) with π = ρ*.)
  **The L = 2 notebook uses the opposite names** ("f" = ρ e_ij, "f̄" = ρ −e_ji); the L = 4 notebook uses Paul's.
- Chain: sites A F (L = 2), A F A F (L = 4). SoV vacuum ⟨0| = e₁ᵀ⊗e₄ᵀ(⊗e₁ᵀ⊗e₄ᵀ). K = Σ e_ij⊗e_ij = |S⟩⟨S| on an A–F pair
  (S = Σ e_i⊗e_i, the invariant), P = permutation. γ := G_c^{−T}.
- Points: h = 1/2, θ = (2/7, −3/5[, 8/11, −11/8]) (numpy checks used θ₃,θ₄ = 0.731, −1.37 in places), z = (3/2, 5/4, 2/3, 4/7).

## References

- Ryan–Volin arXiv:1810.10996 (companion twist frame) covers only the SAME rectangular rep on every site — not f⊗f̄.
- Ryan–Volin arXiv:2002.12341 (any irrep per site): B-operator (2.25), basis (4.12)/(1.2). This is the one that applies.
- Bargheer–Bercini–Lefundes–Ryan arXiv:2603.26617: SoV localisation (G2), free-parameter families W_{κ,ω} (G5)–(G18).

## Results

### L = 2 (`L2-ffbar-SoV-localisation.wb`, 36 cells)

- **qKZ form at the inhomogeneities [exact, symbolic, any twist matrix]:** ℒ_A(0) = hK is rank 1, ℒ_F(0) = −hP, so
  t₁(θ₁) = h ρ₁(Gᵀ)(θ₁₂ − hK), t₁(θ₂) = h ρ₂(G)(θ₁₂ − hK) (ρ_s(·) here = the plain matrix on site s), and
  t₁(θ₁)t₁(θ₂)⁻¹ = ρ(Gᵀ), the global twist rotation.
- **Spectrum [exact]:** simple; 12 of 16 eigenstates are product states (rank 1 as 4×4), 4 zero-weight states mix with the
  singlet.
- **gl(4) B-operator [exact]:** B(u) = Σ_{k₁∈{1,2,3}} Σ_{K₂⊂{1,2,3},|K₂|=2} T_{k₁1}(u) qmin_{K₂;(1,k₁+1)}(u+h)
  qmin_{(123);(1,K₂+1)}(u+2h), untwisted T; degree 12, commuting, simple spectrum, roots = Ryan–Volin GT dual diagonals.
- **SoV basis [exact]:** ⟨a,k| = ⟨0| t_a(θ₁−h) t₁(θ₂)^k, a,k = 0..3 (A site: t_a at θ₁−h; F site: t₁ⁿ at θ₂) — plain transfer
  matrices, no embedding morphism needed. Rank 16, B-eigen, exactly z-independent; wave functions τ_a(θ₁−h)τ₁(θ₂)^k.
- **One-site rotations [exact]:** on any globally rotated vacuum Ω (⟨Ω|K = 0), ⟨Ω|t₁(θ₁)⁻¹ = ⟨Ω|ρ₁(G^{−T})/(hθ₁₂), so a
  one-site twist rotation on top of a global one has overlap hθ₁₂⟨Ω|Ψ⟩/τ₁(θ₁); iterating adds a g-independent ⟨S|-term
  seen only by the 4 mixed states.

### L = 4 (`L2-…wb` cells 30–35 and `L4-ffbar-free-rotations.wb`, 41 cells)

- **qKZ form [exact]:** t₁(θ₁) = h G₁ᵀR₁₄R₁₃R₁₂, t₁(θ₂) = −h R₂₃R₂₄G₂R₂₁, t₁(θ₃) = h R₃₂R₃₁G₃ᵀR₃₄, t₁(θ₄) = −h G₄R₄₁R₄₂R₄₃;
  A–F pairs: R = θ ± hK (scalar off the pair singlet); A–A and F–F pairs: R = θ ± hP.
- **Localisation on a generic globally rotated vacuum Ω [exact]:** end sites only (t₁(θ₁)⁻¹, t₁(θ₄)⁻¹), end blocks
  ⟨Ω|t₁(θ₁)⁻¹t₁(θ₂) = c₁₂⟨Ω|ρ₁(G^{−T})ρ₂(G) with c₁₂ = θ₂₃(θ₂₄−h)/((θ₁₃+h)θ₁₄) (and {3,4}); of the 16 L = 2 SoV words only
  {1, t₃(θ₁−h)} × {1, t₁(θ₂)} localise on block {1,2} (2603-(G2) analogue with one operator per site).
- **On the bare SoV vacuum ⟨0| (g = 1, Paul's choice) [exact]:** EVERY site localises (the singlet term of the interior
  operators vanishes): X = t₁(θ₁)⁻¹ → ρ₁(γ), Y = t₁(θ₂) → ρ₂(γ) (const hθ₁₂θ₂₃(θ₂₄−h)), Z = t₁(θ₃) → ρ₃(γ⁻¹)
  (const hθ₃₂(θ₃₁+h)θ₃₄), W = t₁(θ₄)⁻¹ → ρ₄(γ⁻¹). Site vectors e₁ᵀ(1−κ₁G⁻ᵀ), e₄ᵀ(1−κ₂G), e₁ᵀ(1−κ₃Gᵀ), e₄ᵀ(1−κ₄G⁻¹).
  Overlaps: ∏(1−κ_s ô_s)⟨0|Ψ⟩ for single sites and adjacent pairs; same-type pairs pick up an exchange term from
  R = θ ± hP, e.g. {1,3}: [(1−κ₁ô₁)(1−κ₃ô₃) − (h/θ₁₃)κ₁κ₃(ô₁ô₃−1)]⟨0|Ψ⟩, from ⟨0|ρ₁(γ)ρ₃(γ⁻¹) = ⟨0|[(θ₁₃−h)o₁o₃ + h]/θ₁₃.
  Both ends {1,4} need longer words (singlet-channel obstruction). Exchange dressing on the bare vacuum:
  ⟨0|t₁(θ₁) = hθ₁₂θ₁₄[θ₁₃⟨0|ρ₁(Gᵀ) + h⟨0|ρ₃(Gᵀ)].
- **Global torus [exact]:** ρ(X), X ∈ C[G_cᵀ], commutes with all t_a (character ∏f(z_i)^{μ_i}). Parameter count of κ's +
  torus: Jacobian rank 7.
- **No higher twist powers [numeric scan]:** of 2288 fused T_λ (|λ| ≤ 6) at θ_α + mh, none localises as G^{±2}, G^{±3} on
  a site; each site moves only in span{vac, vac·γ^{±1}}.
- **Reach [exact]:** global rotations see only the 84; the 16 κ-products see all six irrep types; κ + torus span 193/256.

### Common rotation of a block — {1,3} [numpy only, 10⁻¹⁴; `checks/chk_common13.py`] — **newest, not yet in a notebook**

Paul's expectation is right: with t₁(θ₁) (not its inverse) the site powers match. On ⟨0| (R_s = ⟨0|ρ_s(Gᵀ),
R₁₃ = ⟨0|ρ₁(Gᵀ)ρ₃(Gᵀ)): ⟨0|t₁(θ₃) = l₃R₃, ⟨0|t₁(θ₁) = α(θ₁₃R₁ + hR₃) with α = hθ₁₂θ₁₄, and ⟨0|t₁(θ₁)t₁(θ₃) = m·R₁₃ exactly
(pure two-site rotation; m = 0.056461 at the point, closed form not fitted). Hence

  ⟨0|(1 + a t₁(θ₁) + b t₁(θ₃) + c t₁(θ₁)t₁(θ₃)) = ⟨0|ρ₁(g)ρ₃(g),  g = 1 − κG_cᵀ,
  a = −κ/(αθ₁₃), b = −κ(θ₁₃−h)/(θ₁₃ l₃), c = κ²/m,

a common one-parameter rotation of sites 1 and 3 (2603 resolvent structure), overlap (1 + aτ₁(θ₁) + bτ₁(θ₃) + c τ₁(θ₁)τ₁(θ₃))⟨0|Ψ⟩.
The same span also gives independent (1−κ₁Gᵀ) ⊗ (1−κ₃Gᵀ). Earlier "no common rotation" statements (built from t₁(θ₁)⁻¹,
which gives G⁻ᵀ on site 1 against Gᵀ on site 3) are superseded.

### Globally rotated reference states [numpy only, 10⁻¹⁴; `checks/chk_rotref.py`]

- The {1,3} common rotation is a global rotation: ⟨0|ρ₁(g)ρ₃(g) only sees the first row of g, so g → g′ = 1 − κ(c₁E₁₁ + E₁₂)
  ∝ exp(−κ′E₁₂), which fixes the F vacuum e₄ᵀ; hence ⟨0|ρ₁(g)ρ₃(g) = ⟨0|ρ(g′) on all four sites (in the 84). (Not with g itself:
  G_cᵀ contains E₃₄, which moves e₄ᵀ.) The {2,4} analogue is generated by E₃₄.
- **Criterion:** a global rotation v⊗w⊗v⊗w (v·w = 0) is as good a reference as ⟨0| — all four sites localise, both common
  rotations exist, same coefficients a, b, c — iff w G vᵀ = 0, i.e. the conjugated twist keeps (G)₄₁ = 0. Checked on random v, w.
- **Reach:** with v = e₁ᵀp(Gᵀ), w = e₄ᵀq(G) and (Gᵏ)₄₁ = 0 for k ≤ 2, (G³)₄₁ = 1: global iff deg p + deg q ≤ 2, good reference iff
  ≤ 1. So ⟨0|P₁₃(κ₁)P₁₃(κ₂) = (e₁ᵀ(1−κ₁Gᵀ)(1−κ₂Gᵀ) ⊗ e₄)^⊗2 (Gᵀ² on a pair of sites) and ⟨0|P₁₃(κ₁)P₂₄(κ₂) = one global
  rotation on all four sites (A: e₁ᵀ(1−κ₁Gᵀ), F: e₄ᵀ(1−κ₂G)), overlaps p₁₃(τ)p₂₄(τ)⟨0|Ψ⟩ — polynomial, no resolvent. Degree-2
  states are not good references; a third step leaves the global orbit.

### Rotations of the 2-parameter state [numpy; `checks/chk_span13_words.py`, `checks/scan_*.py`]

- With X = t₁(θ₁), Y = t₁(θ₃): the nine ⟨0|X^aY^b (a, b ≤ 2) are independent and span exactly
  {e_i⊗e₄⊗e_j⊗e₄ : i, j ≤ 3}. The global rotations inside are v⊗e₄⊗v⊗e₄, v ∈ span(e₁,e₂,e₃) — exactly
  ⟨0|P₁₃(κ₁)P₁₃(κ₂), generated by E₁₂, E₁₃ (+E₁₁) of the gl(3) on {1,2,3}:
  e₁ᵀ(1−κ₁Gᵀ)(1−κ₂Gᵀ) = (1 − sc₁ + p(c₁²−c₂))e₁ + (pc₁ − s)e₂ + pe₃ (s = κ₁+κ₂, p = κ₁κ₂). (Gᵀ)² on site 1 alone is also in
  the span (combinations of words beat the single-T_λ "no higher powers" scan); v₄ ≠ 0 is not.
- The 2-parameter state Ω₂ fails w G vᵀ = 0 and behaves exactly like a generic global rotation: only end-growing words
  localise, as fixed twist rotations — t₁(θ₁)⁻¹ → ρ₁(γ), ·t₁(θ₂) → ρ₁₂(γ), ·t₁(θ₃)⁻¹ → ρ₁₂₃(γ), ·t₁(θ₄) → ρ(γ) (torus);
  mirror from site 4 with γ⁻¹. Free parameters only on one end site. Table in `ABJM-3pt-transfer-matrices.md` §4.

### ABJM three-point functions (diagonal frame) — **see `ABJM-3pt-transfer-matrices.md`**

Combines the ABJM note (`StructuringABJM.pdf`), arXiv:2603.26617 and the rotation results. Frame map
t^{G_c}R = R t^Λ (R = V^{−T} on A, V on F). The note's (3.6) = ⟨W₋₁^L| with the first pairs rotated by 1 − Λ/z₄
(W_a = (e₁+ae₄)⊗(e₄−ae₁)). At L = 2: odd bridge (one field) works exactly, C ∝ (1 − z₄/(c₁τ₁(θ₁)))/(1 − z₄/z₁)·⟨W₋₁^L|Ψ⟩, checked on on-shell HW
states through the untwisting limit; explicit C_ℓ table for all L = 2 HW states (adjoints seen only by odd bridges,
singlets only by even); even bridge
(a pair) is not reached by localising words — the missing piece is Λ on the F site alone; the highest-weight trick
(free global e^{bE₁₄}) gives an exact finite-twist formula with b = z₁/(z₁−z₄), but b diverges and the untwisting limit
is wrong. 2603 needs (and has) a continuous κ-family for this; for ABJM pairs it is the open problem.

## Open / next

1. **ABJM pair κ-family** (main blocker, `ABJM-3pt-transfer-matrices.md` §6): an F-site-only twist rotation at the left
   end, or a word combination cancelling the singlet in ⟨W|t₁(θ₂); try L = 3 pairs. Then K_κ for ℓ pairs and the
   untwisting limit with κ → 1.
2. Odd bridges: homogeneous limit of the one-field formula; ℓ = 3 fields. Integrability of ⟨W_a^L| (note's (4.8)).
3. Put the {1,3} common rotation into a notebook exactly (closed form of m; second (h,θ) point).
4. More parameters per site: only one twist power per site localises from single T_λ; combinations reach (Gᵀ)² on a
   site (3⊗3 span) but not degree 3.
5. L = 4 SoV basis; Q-function (Slater determinant) form of the overlaps and of ⟨0|Ψ⟩ / ⟨W^L|Ψ⟩.

## Files

| File | Content |
|---|---|
| `L2-ffbar-SoV-localisation.wb` + `-spec.md` | L = 2: qKZ, spectrum, B-operator, SoV basis, one-site rotations; L = 4 preview (cells 30–35). Summary cell 36. |
| `L4-ffbar-free-rotations.wb` + `-spec.md` | L = 4: localisation inventory, free-parameter families on ⟨0|, torus, fused-T scan, reach. Summary in the last cells. |
| `dual-lax-t3.wb` | t₂, t₃ as single-pass transfer matrices. t₃: fused ∧³ Lax from quantum minors (symbolic), t₃ = ∏g_s·tr[ℒ̄₁⋯ℒ̄_L det G G^{−T}] (exact, L = 2, 4, four twists), why the naive dual Lax fails (relative A/F shift 4h), vector form, L = 6 commutation test. t₂ (§4): fused ∧² Lax (6-dim aux, A: g = u−θ−h, shift 0; F: g = u−θ, shift +h), t₂ = ∏g_s·tr[ℒ⁽²⁾₁⋯ℒ⁽²⁾_L ∧²G] (exact, L = 2, 4, four twists), vector form. |
| `abjm-chain.wb` | The ABJM chain: definition, Minahan–Zarembo Hamiltonian (arXiv:0806.3951, Bak–Rey 0807.2063) with sanity checks, exact [H_MZ, t̄_a(u)] = 0 (a = 1, 2, 3, symbolic u, L = 2 pairs) at θ_A − θ_F = 2h, uniqueness of that point (symbolic θ_A − θ_F), H_MZ = U(0)⁻¹U′(0) + const from U = t̄₁(θ_F+u) t̄₃(θ_A−h+u), numerical check at L = 3 pairs. |
| `abjm-LR-symmetry.wb` | 2603's LR symmetry for ABJM: Z₂'s 𝒫 (t̄₁(u) ↔ t̄₃(u+h)) and 𝒬 (t̄₁(u) → t̄₁(2h−u)); all HW cyclic t-eigenstates at L = 2, 3, 4 pairs (6, 26, 152), C_ℓ for all bridges, vanishing rules (R-charge, bridge parity for L+p odd, no LR selection rule, isolated zeros). States saved in `data/abjm-hw-states.wxf`. |
| `abjm-Qfunctions.wb` | Q₁, Q₁₂, Q₁₂₃ of all saved states (Baxter equation in gauge d = V₄, Casoratians / vacuum Casoratians; degrees = magnon numbers for all 184 states). LR in Q-language: Q₁(u) ∝ Q₁₂₃(u+h); agrees with the τ criterion for every state; every non-LR state has an LR partner. Saved in `data/abjm-Qfunctions.wxf` (+ `data/abjm-tau-L4.wxf`). |
| `abjm-3pt-configs.wb` | Search for a 3-pt configuration with an LR selection rule: all two-block BPS bras (6 GL(4) orbits × bridge length × A/F start, L = 2, 3, 4): none at L = 3, 4. |
| `abjm-LR-factorized.wb` | Does the LR-annihilated space (left kernel of t̄₁(u) − t̄₃(u+h)) contain factorised bras? L = 2, untwisted and twisted (z₁z₄ = z₂z₃ = 1): only trivial eigen-bras. |
| `abjm-other-symmetries.wb` | Parity-type symmetries: LR, P (τ̄₁(u) = τ̄₁(2h−u) ⇔ Q₁ ∝ Q₁₂₃(2h−u)), LP (τ̄₁(u) = τ̄₃(3h−u) ⇔ Q₁ ∝ Q₁(h−u)) on all saved states; product states preserving them (P: all two-site states). |
| `data/abjm-tau.wxf` | τ̄₁, τ̄₂, τ̄₃ eigenvalue polynomials of all saved states (L = 2, 3, 4). |
| `data/abjm-hw-states.wxf` | Saved verified states (load with `Import[…, "WXF"]`; key "description" documents conventions). |
| `abjm-bridges.wb` | Port of `abjm_bridges*.py`. §1–2 at θ = 0 (superseded: singular point, S not in the t-algebra, twisted-shift workaround). Shift test U = t₁(θ_F)t̄₃(θ_A−h) ∝ S; Minahan–Zarembo H commutes with the t's only at **θ_A − θ_F = 2h (ABJM point)**. §3: on-shell states, C_ℓ table and odd-bridge untwisting test at the ABJM point. Set `Lp` for the number of pairs. Summary in the last cell. |
| `ABJM-3pt-transfer-matrices.md` | ABJM 3-pt functions from localising transfer matrices: sources (note, 2603), dictionary, frame map, L = 2 results, untwisting test, open problems, script index. |
| `checks/*.py` | numpy cross-checks (run with `python checks/<name>.py`; scripts `exec` siblings in the same folder): `ffbar.py`–`ffbar4.py` qKZ fits and end-site localisation; `check5.py` end blocks; `check6.py` SoV basis z-independence; `g2test.py`, `g2inv.py` (G2)-type localisation inventory; `chk7.py` interior sites at g = 1; `chk13.py` {1,3} exchange identity; `chk_common13.py` {1,3} common rotation; `chk_rotref.py` rotated reference states, criterion, iterated common rotations (`rotbase.py` = its shared setup); `chk_span13_words.py` 3⊗3 word span; `scan_*.py` Ω₂ localisation scan; `chk_frame.py` companion↔diagonal frame map; `abjm_*.py` ABJM 3-pt localisation, untwisting and explicit bridge overlaps (index in `ABJM-3pt-transfer-matrices.md` §7). |
| `ref/*.txt` | cell-source dumps of two gaudin-cg notebooks the specs cite (`L2-gl4-rect-untwist.wb`: gl(4) companion matrix, P_aux, quantum minors; `L2-rotation-sov.wb`: gl(3) L = 2 B-operator and SoV basis). |

Notebook caveats: both .wb files were built as JSON and run in a kernel by an agent — **no native outputs**; results are in the
markdown cell under each code cell. Rerun in VS Code to store outputs. The L4 notebook's four `Export`s write to
`$TemporaryDirectory` (`ffbar-*.m`). The specs mention a "session scratchpad": its scripts are now `checks/`, its cell dumps
`ref/`; the paper text extracts were not copied (use the arXiv numbers above).
