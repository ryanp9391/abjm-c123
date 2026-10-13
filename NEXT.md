# Next — abjm-c123

*Updated 2026-10-13 (morning)*

**Focus:** tree-level structure constants in ABJM theory, full su(4) sector. Active thread: the su(4) alternating
spin chain (fund ⊗ anti-fund ⊗ …) — SoV basis, localisation, free-parameter rotations on few sites — now being
connected to the ABJM 3-pt function of the note `StructuringABJM.pdf` (repo root) following arXiv:2603.26617.
Everything is in `Paul/Mathematica/su4-alternating/`: **read its `README.md`, then `ABJM-3pt-transfer-matrices.md`**.

## RESUME HERE (session of 2026-10-12, ended with the Wolfram licence expired)

**Complete notes (2026-10-13): `Paul/Mathematica/su4-alternating/notes/complete-notes/abjm-complete-notes.pdf` (61 pp.) — read this
first; it supersedes the 12-Oct tex note and condenses everything below. Open problems in its §19.**
**Middle even bridges (2 ≤ M ≤ L−2, first L = 4, M = 2) are still OPEN** — the κ-family covers M ∈ {1, L−1, L}. New untested
route (notes §19 item 1): Proposition 1 works pair by pair, so the companion bra with the first M pairs of (e₁⊗ũ)^L rotated by
1−κG_c (A: e₁ᵀ(1−κG_cᵀ), F: ũ(1−κG_c)⁻¹) has the overlaps of ⟨W_{a′}^M⊗W_{κ₀}^{L−M}| on the sectors — no κ-family, no z₂ = z₃;
needs the closed form of its SoV coefficients (Proposition-2 analogue). Test at L = 2, M = 1 first, then L = 4, M = 2.

0. **Twisted-vacuum overlap — RESOLVED 2026-10-13.** ⟨W_a^L| is not itself a twisted vacuum, but ⟨0_w|Ψ⟩ = (v₁u₄)^L⟨W_a^L|Ψ⟩ for every Ψ of weight (μ,0,0,−μ), with w = (w₁, z₂, z₃, w₄) and a = −(z₁−z₂)(z₁−z₃)(w₁−z₁)/((z₄−z₂)(z₄−z₃)(w₁−z₄)); ⟨W_a^L|Ψ⟩ = 0 on all other weights
   [exact, L = 2]. ⟨0_w|Ψ⟩ is a finite SoV Q-sum 𝕎_w (agent B, 1e-44), so the even bridge M = 1 is in Q's at finite twist:
   ⟨W_{a′}⊗W_a^{L−1}|Ψ⟩ = P_κ(τ)(v₁u₄)^{−L}⟨0_w|0̄⟩⟨0|Ψ⟩𝕎_w(Ψ)/(1−κz₁)². Read
   **`Paul/Mathematica/su4-alternating/notes/w-a-reference-status.md`** first (frames, reference bras, agents A/B/C, §5 formula).
   Agents of 2026-10-13 (one shared kernel; it faulted once from agent C's runaway memoised cell and was restarted 08:03):
   A (`reference-admissibility-Kkappa.wb`): changing the reference to ⟨0_w| fails (no κ-family at z₂ ≠ z₃), superseded by the
   overlap identity. B (`Wa-sov-expansion-L4.wb`): ⟨W_a^L| dense in SoV words, no determinant; 𝕎_w Q-sum verified.
   C (`even-bridge-epsilon-expansion.wb`, ported: exact −1/√40, 1/√24 at L = 2, 1e-14 at L = 3; ε-scans python): C_{M=1} = K₁(τ_Ψ₀)·Ω(Ψ₀), K₁ = 2 − A − B − C untwisted Q's,
   Ω bridge independent (= ε² coefficient of the 𝕎_w overlap); no even/odd ratio relation.
   **B follow-up DONE:** w₁ → ∞ gives ⟨W_a^L|Ψ⟩ = (a/κ₀)ⁿ(λ₁u₄)^{−L}⟨0|∏_F ρ_F((G_c−z₂)(G_c−z₃))|Ψ⟩ = F-only SoV sum,
   at L = 2 a bidegree-(2,2) polynomial in τ₁(θ₂), τ₁(θ₄) [1e-51]; the L = 2 even bridge from it matches 1/√40, 1/√24 within
   1e-11 at the ABJM point. Q-complete up to the SoV norm ⟨0|Ψ⟩/‖Ψ‖. Status note §5.
   **Agent D DONE** (`Wa-Fonly-L6.wb`): F-only form exact at L = 3 (27 words); closed form FACTORISES,
   ⟨W_a^L|Ψ⟩/⟨0|Ψ⟩ = (a/κ₀)ⁿ(λ₁u₄)^{−L}B_{z₂}(y)B_{z₃}(y), B_z = Σ_S(−z)^{L−|S|}∏_{s<t∈S}ν_st²/(ν_st²−h²)∏_{s∈S}y_s, y_s = τ₁(θ_s)/κ_s
   [exact L = 2, 3; model L = 4]; L = 3 even bridge |C| = data to 8 digits on the six [1,0,1] states (norm from eigenvector).
   **Agent E DONE** (`abjm-sov-norm-L4.wb`): the L = 2 even bridge is FULLY Q-based, |C| = lim 𝒩_pre·𝔸·𝕎/√𝔹 with
   𝔹 = Σ τ̃Mτ from regularised bra/ket eigenvalue data; (𝔸,𝕎,𝔹) → (10,25,2.5e6), (2,9,7776) = 1/√40, 1/√24 (4e-16).
   **Agent F DONE** (`even-bridge-Qform.wb`, consolidated, 58 cells): all ingredients re-verified; L = 2 analytic limit
   |C| = K₁·B₁²/√𝔹₀ with K₁ = 2 − 8τ − 1/(2τ) + 64xτ, B₁ = 1 + 4τ + 2τ′ + 16τ² (τ = τ₁(θ_F), τ′ = ∂_uτ₁(θ_F), x = τ₁(θ_F+h));
   𝔹 = the SoV norm is KNOWN to Paul — never derive it. **Next:** L ≥ 3 untwisted limit (only the product of factors is finite); proof of the B_z form; middle bridges (L = 4, M = 2); rerun
   the unrun notebooks of 2026-10-12 one at a time — STOPPED by Paul ~11:50 after 4 of 10 (log
   `Paul/Mathematica/su4-alternating/notes/rerun-2026-10-13.md`): done pair-group-rotation (cell 17: four 2-factor words localise with
   constants 51/31, 71/31 × ∏λ_k, not ∏λ_k; K_κ unaffected; L = 3 rows skipped as too heavy), even-bridge-lpairs (F1/F2 labels swapped
   in yesterday's table, same limits), odd-bridge-general-L (all reproduced), sov-Q-dictionary-L4 (saved, not logged).
   sov-abjm-point-limit.wb was mid-run (cells 27+ cancelled; earlier cells may hold unsaved outputs). Not started:
   even-bridge-sov-formula, sov-onepair-determinant, sov-twoblock-expansion-L4, twoblock-polarisation-forms,
   twisted-vacuum-localisation-L2, twisted-vacuum-companion-check.
1. **Read `Paul/Mathematica/su4-alternating/notes/session-2026-10-12-agent-reports.md`** (all results of the day, condensed)
   and items 2h–2m below. Headline: the pair κ-family / even bridge exists (2k: exact 5-word K_κ at the ABJM point with
   z₂ = z₃); SoV at the ABJM point is regularised (2m); SoV dictionary + dual basis verified (2l).
2. **Renew the Wolfram Engine licence**, then rerun today's notebooks in VS Code (all have unrun cells with python
   results in markdown; small port bugs possible). Python fallbacks: `checks/session-2026-10-12/`, `checks/abjmlimit/`,
   `checks/twoblock-polforms/` (each with RESULTS.md or a log).
3. All agents finished. `even-bridge-sov-formula.wb` (notes §9d): K_κ is NOT a common rotation on ⟨0|_comp and the
   singlets are invisible from ⟨0| with pair-localised words (F site never gets e₁); general SoV-normalised formula
   |C_K| = |P_K(τ_Ψ)/P_K(τ_vac)|·|λ^{w−w_vac}𝒩(Ψ)/𝒩(vac)|^{1/2} verified at L = 2; regularised ket basis and Q-forms of x, w, y
   at the ABJM point; caveat: z₂ = z₃ has no SoV basis (gl(2) symmetry) — combine with z₂−z₃ = O(ε). Done: `even-bridge-lpairs.wb`
   (notes file §9c): no common-κ ℓ-pair word — the structure is pretwisted (end-growing); the one-pair word gives all
   even bridges with M ∈ {1, L−1, L} (C_M = (−1)^{L+μ₄}C_{L−M}), verified O(ε) on all 26 L = 3 states; z₂−z₃ = O(ε)
   suffices; per-pair factor explicit in Q's (τ̂₃(u₀) = dual-Q ratio). **First open case: L = 4, M = 2** (needs a second
   integrable family or a different reference). Done: `pair-rotation-locality-L3.wb` (generic θ, z: even bridge NOT locally reachable —
   consistent with 2k, the ABJM point + z₂ = z₃ is essential; notes file §9b).
4. Next: (i) the middle bridges 2 ≤ M ≤ L−2 (L = 4, M = 2 first); (ii) ⟨W_a^L|Ψ⟩ in SoV/Q form (2603's 𝕎); (iii) the norm
   / measure of the mixed chain (2603's 𝔹); (iv) relax z₂ = z₃; (v) write-up (update the tex note).
5. Nothing of today is committed (git status: 20+ new notebooks, notes, checks, data/abjm-twoblock-forms.json).

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
   **Complete write-up of the whole thread** (conventions + RV dictionary, t₂/t₃, ABJM point, localisation, 3-pt
   construction, states/Q's, symmetries, selection rules, SoV, open problems, notebook/data index, superseded
   statements): `su4-alternating/notes/abjm-transfer-matrices-symmetries.tex/.pdf` (9 pp., 2026-10-12).
2h. **Twisted vacuum as reference (2026-10-12, `twisted-vacuum-*.wb`, five notebooks, Sonnet agents):** companion-frame
   eigenbra of t^{G_c(w)} with the vacuum eigenvalue is v_w ⊗ u_w ⊗ …, v_w = (w₁³,w₁²,w₁,1), u_w = row 4 of V(w)⁻¹.
   Acting with t^z: obstruction q = u_w G_c(z) v_w = −∏(w₁−zᵢ)/∏_{m≤3}(w₄−w_m). Generic w: only end sites (inverse
   words, fixed G^{∓T}) and end-growing chains (fixed, through all 6 letters at L = 3 pairs); the only free families on
   pair 1 are v(α+βG⁻ᵀ)⊗u and vG⁻ᵀ⊗u(α+βG). Locus w₁ ∈ {z}: A sites frozen (G_z-eigenvectors), F sites free: common
   (1−κG_z) on the F sites of the first ℓ pairs from words in t₁(θ_{2k}) (closed-form coefficients, L = 3 pairs).
   Diagonal frame + untwisting (κ → 1 with z → 1): rotated pairs → e₁⊗(e₄+αe₂+βe₃) (locus) — both O₂–O₃ pairings
   vanish, overlaps only the BPS vacuum; generic w with w₂,w₃ = z₂,z₃ → e₁⊗(e₄+α̃e₁), sees [L−1,0,L−1] only. Neither
   gives the note's W₋₁ (needs an O(1) A-site rotation e₁ → e₁−e₄; all A rotations available are fixed twist powers).
   Notebooks have no stored outputs (shared kernel) — rerun in VS Code.
2i. **Why 2603 escapes this and what it means for ABJM (2026-10-12, second half):** 2603's site is the 6 = ∧²4, vacuum
   Z = e₁∧e₂; its κ-family is the one-point word 1 − κt₊(θ₁) + κ²t₀(θ₁) (t₀ = ∧² transfer matrix) = the GROUP element
   ∧²(1−κg) on the site; localisation (their App. K) needs no goodness condition, the rotated vector stays null and its
   pairing with the unrotated one is quadratic in D = 1−κg, so non-extremal. ABJM analogue needed: a COMMON pair
   rotation (A by D = 1−κΛ, F by D⁻¹ ∝ adj D, cubic). It is automatically BPS, its pairings with the reference are
   −κ u·Λv and Σvⱼuⱼ/(1−κzⱼ) (O(1) for non-good references), and on W_a^L it gives pair 1 → W_{a′} with a′ free, so a
   finite HW shift lands on vac^{L−1}⊗W₋₁ (the note's bra) with a continuous limit. Status: `pair-rotation-wordsearch.wb`
   (2 pairs, W_a², all t₁,t₂,t₃ shifts ≤ 2 factors): the localised-on-pair-1 span contains the whole e₁/e₄ block, so
   T_A is reachable but only with ≥ 3 words (t₁ at θ₁−2h, θ₁−3h + mixed), non-uniquely — at 2 pairs this is nearly
   trivial (commuting algebra is rich; generic reference spans all 256). **Resume here:** locality test at L = 3 pairs
   (does a few-word combination with universal coefficients give W_{a′}⊗W_a⊗W_a from W_a³?), and a 2603-App.-K-style
   derivation for the pair at u = θ_A / θ_F, where the pair vacuum has the same 6 + 2 "bonus" zero pattern
   (⟨e₁e₄|(L_AL_F)_{jk} = 0 for j > k, and (j,4), j = 2,3 at u = θ_A; (1,k), k = 2,3 at u = θ_F) — ideally at the ABJM
   point θ_A − θ_F = 2h.
   **Done exactly (`abjm-nongood-reference-limit.wb`, 2 pairs):** F-field one-field bridge formula
   ⟨W_a^L|t₁(θ₁)⁻¹(α + βt₁(θ₂)) = ⟨W_a^L|ρ₁(Λ⁻¹)ρ₂(αc₁ + βc₁₂Λ), c₁₂ = θ₂₃(θ₂₄−h)/((θ₁₃+h)θ₁₄); with α+β ∝ ε the
   limit bra is explicit (pair 1 → (e₁+ae₄)⊗(e₄−a′e₁), a′ = a(c+φ₁)/(c+φ₄)), finite HW shift gives vac^L with one F
   site = e₁, the 0·∞ product converges like ε, φ- and c-independent, to (a−a′)C₁. §5.5 is the κ = z₄ point of the
   A-site family. Even bridge not reached by these families or the end-growing chain (confirmed).
2j. **Odd bridges at general L — closed forms (2026-10-12, `odd-bridge-general-L.wb`, proved in §1, exact L = 2, 3,
   1e-16 at L = 4):** on any reference (v⊗u)^L with v·u = 0 (same pair everywhere), the end-growing chain
   P_{k−1}t₁(θ_k)^{∓1} = d_k P_k (A sites inverse letters → Λ⁻¹, F sites forward → Λ) with
   1/d_k = h∏_{s≠k}ϑ_ks (k odd; ϑ_ks = θ_ks for F, θ_ks+h for A), d_k = −h∏_{s≠k}ϑ̃_ks (k even; ϑ̃_ks = θ_ks for A,
   θ_ks−h for F), c_k = d₁⋯d_k; c₁ = §5.5's formula (now proved), c₁₂ = c₂, c_{2L} = ∏_{F<F'}(h²−θ²_FF')/∏_{A<A'}(h²−θ²_AA').
   Both one-field families (A-site 1−κΛ⁻¹ via (α+βt₁(θ_k)⁻¹), F-site via (α+βt₁(θ_k)) after the chain) localise
   exactly at every placement k = 1..2L (L = 3) and untwist (rate O(ε), φ/c-independent) to V₀ + (a−a′)C₁: the ℓ = 1
   odd bridge sees only [L−1,0,L−1]. Homogeneous ABJM-point limit is REGULAR (c₁ = (2h²)^{−L}, c₁₂ = 1); the singular
   locus is θ_A = θ_F; normalised τ̂_k = d_k^{±1}τ(θ_k) stay finite, formula
   ⟨bra_ε|Ψ_ε⟩ = [∏_{j<k}τ̂_j^{(−1)^j}](cε−1+τ̂_k)⟨W_a^L|Ψ_ε⟩. Notebook cells 16–34 unrun (kernel died).
   **BLOCKER (2026-10-12 evening): the Wolfram Engine licence has expired** (wolfram.exe: "license has expired, no
   valid password found"); the shared kernel died after the VS Code extension host reached 14–19 GB. All notebooks from
   the second half of the day were written but not (fully) run; results are in markdown cells / python in the
   session scratchpad. Renew the licence, then rerun.
2k. **BREAKTHROUGH — the pair κ-family (even bridge) exists at the ABJM point (2026-10-12 night, `pair-group-rotation.wb`,
   python in `checks/session-2026-10-12/pairgroup/`):** on ⟨W_a^L| with pair-periodic θ, θ_A − θ_F = 2h, and a twist with
   z₂ = z₃ =: y, the common group rotation of pair 1 is the exact 5-word combination (X₁ = t₁(θ₂), X₂ = t₂(θ₂+h),
   X₃ = t̂₃(θ₂+h) reduced, s = 1−κ(z₁+z₄), λ_k L-dependent products of θ-differences, see pgbase.py)
   K_κ = s² + κs[X₁/(h²λ₁) + X₂/(2h⁴yλ₂) + X₃/(6h⁵y²λ₃)] + κ²X₃X₁/(24h⁷y²λ₁λ₃),
   ⟨W_a^L|K_κ = (1−κz₁)²⟨W_{a′}⊗W_a^{L−1}|, a′ = a(1−κz₄)/(1−κz₁). Unique, independent of a, local (exact at L = 2, 3).
   Mechanism = 2603 App. K: the ∧^k structure lives on the fused pair; t₂(θ_F+h) localises ONLY at the ABJM point
   (hence all generic-θ searches failed); z₂ ≠ z₃ leaves a zero-weight remainder ∝ z₃E₂₂+z₂E₃₃ (6 directions vs 5
   operators). Untwisting: κ = 1/z₄ gives ⟨vac⊗W_a^{L−1}|Ψ⟩ = (z₁²−z₁Σ+Π)/(z₁−z₄)²⟨W_a^L|Ψ⟩ → data C (L = 2 singlets
   −1/√40, 1/√24); κ → 1 with a′−a = −1 → the note's vac^{L−1}⊗W₋₁ at L = 3 (O(ε)). Fails on generic BPS references.
   **In progress:** ℓ-pair generalisation + all M (`even-bridge-lpairs.wb`), the same word on the companion SoV
   vacuum ⟨0| to get C_M = P_κ(τ)·⟨0|Ψ⟩/norm in Q-form (`even-bridge-sov-formula.wb`).
2l. **SoV results (same night):** dictionary verified at 4 sites, no renormalisation needed (`sov-Q-dictionary-L4.wb`,
   `notes/rv-sov-dictionary.md`): F sites single-Q ratios (any j), A sites 3×3 Casoratians with one column omitted
   (any triple); dual ket basis |0̄⟩ = e₄e₁e₄e₁ with A: t₁(θ_A)^m, F: t_n(θ_F+nh); Gram matrix triangular (charge
   conservation), no diagonal measure; 1/𝒩 = Σ τ̃_y M_{yx} τ_x. One pair: ⟨v⊗u|Ψ⟩ = ⟨0|Ψ⟩det[w; Q₁;Q₂;Q₃]/Q₁₂₃ exactly
   (`sov-onepair-determinant.wb`); two pairs never factorise; the F–F exchange terms are NOT a separable (Vandermonde)
   measure — no one-row-per-site determinant exists in the τ-language; a 2603-type determinant needs the (unknown)
   right SoV basis / measure of the mixed 4/4̄ chain. ⟨0|_comp = λ₁²vac^Lρ(n′D′) is in the HW orbit of the
   physical vacuum (`sov-twoblock-expansion-L4.wb`): the M = 1 even bridge at L = 2 is a 36-word SoV sum = bordered
   4×4 Casoratian, untwisting matches data; the SoV basis is block-diagonal in the charge Σ(i−1)_A + Σ(4−j)_F.
   Polarisation forms (`twoblock-polarisation-forms.wb`, `data/abjm-twoblock-forms.json`): F is a polynomial in 8
   Borel invariants with exponents fixed by [p,q,r]; the note's W₋₁ polarisation is BLIND to all states with L+p odd
   (F/F_vac ∝ (x+y), x = −y forced) — the even-bridge parity rule is a property of the polarisation — and to all q > 0
   states (need the minors mv, mu). The [2,0,2] E = 16 state at L = 3 vanishes for every two-block polarisation.
2m. **SoV at the ABJM point — solved (`sov-abjm-point-limit.wb`, `checks/abjmlimit/`):** at u₀ = θ_A − h = θ_F + h the exact
   identities t₃(u₀) = 0 and t₂(u₀) = t₁(u₀)t₁(θ_F) hold (per resonant pair the generic words drop to 9: rank 81 at 4
   sites with pair-by-pair θ); B has Jordan blocks (10 simple + 3 J₂ per pair). Regularised basis = d → 2h limit of the
   RV basis: per resonant A site O(a,k) = {1, t₁(u₀), W = t₂′(u₀) − t₁′(u₀)t₁(θ_F) (k ≤ 2; t₂(u₀) at k = 3),
   Z = t₃′(u₀) = g₃′t̂₃(u₀)}, rank 256, wave functions ∏_pairs λ(a,k)τ_F^k with λ ∈ {1, x = τ₁(u₀), w, y = τ₃′(u₀)};
   the missing A-site data is in the u-derivative of the twisted Baxter equation. Fixed-bra coefficients diverge like
   δ^{−(#resonant pairs)} but the overlap converges (verified 10⁻⁴⁵). Caveat: FULLY pair-periodic θ (coincident
   inhomogeneities) collapses the generic words to rank 70 and needs derivatives to order 7 — the homogeneous limit is a
   separate open problem; use pair-by-pair θ = (θ_F1+2h, θ_F1, θ_F2+2h, θ_F2, …). The K_κ word of 2k is built from
   exactly these regularised pair words (X₂ = t₁(u₀)X₁, X₃ ∝ Z).
3. Integrability of ⟨W_a^L| (note's (4.8)); later the ⟨W^L|Ψ⟩ determinant (Q stage).
4. Notebook versions of the numpy results ({1,3} common rotation closed form, frame map, odd-bridge formula).

## Housekeeping

- The two `.wb` notebooks have no native outputs — rerun in VS Code to store them.
- Committed up to 56b9407 (thread + `chk_rotref.py`); the ABJM write-up, README/NEXT updates and the new `checks/`
  scripts (`rotbase.py`, `chk_span13_words.py`, `scan_*.py`, `chk_frame.py`, `abjm_*.py`) are **not committed yet**
  (all smoke-tested, exit 0).
- 2603 LaTeX source is not in the repo (read from arXiv e-print 2603.26617v1, `short.tex`).
