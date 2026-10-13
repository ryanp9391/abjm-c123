# Session 2026-10-12 — consolidated agent reports (twisted vacuum → pair κ-family → SoV)

Verbatim-condensed final reports of the Sonnet/Opus agents of this session, kept here because the Wolfram Engine
licence expired mid-session (kernel dead): most notebooks are written but UNRUN (python results pasted as markdown
under each code cell); the python that produced the numbers is in `checks/session-2026-10-12/`, `checks/abjmlimit/`,
`checks/twoblock-polforms/`. Evidence tags: [exact] rational/symbolic, [hp] mpmath/CRT high precision, [num] double.
Conventions as in README (A F A F, bras row vectors, A: v → v g, F: u → u g^{−T}, Lax L_A = (u−θ)+hK, L_F = (u−θ)−hP,
T = L₁⋯L_{2L} site 1 leftmost, twist on the right, companion G_c(z), diagonal Λ, R = ⊗(V^{−T} on A, V on F)).
Standard point h = 1/2, θ = (2/7, −3/5, 8/11, −11/8[, 5/13, −7/9]), z = (3/2, 5/4, 2/3, 4/7), w = (7/5, 9/8, 3/7, 5/9).

## 1. Twisted vacuum as reference (`twisted-vacuum-localisation-L2.wb`, `-L4.wb`, `twisted-vacuum-blocks-L3.wb`)

- ⟨0_w| := eigenbra of t^{G_c(w)} with the vacuum eigenvalue = v_w ⊗ u_w ⊗ … with v_w = (w₁³,w₁²,w₁,1) (right eigenvector of
  G_c(w), eigenvalue w₁), u_w = row 4 of V(w)⁻¹, (u_w)_k = (−1)^{k−1}e_{k−1}(w₁,w₂,w₃)/∏_{m≤3}(w₄−w_m) (left eigenvector,
  eigenvalue w₄; depends on w₄ only through normalisation); u_w·v_w = 0. Exact eigenbra of t₁,t₂,t₃ [exact]; eigenvalue at
  4 sites: w₁(u−θ₁+h)(u−θ₂)(u−θ₃+h)(u−θ₄) + (w₂+w₃)∏(u−θ_s) + w₄(u−θ₁)(u−θ₂−h)(u−θ₃)(u−θ₄−h).
- Obstruction q(w,z) = u_w G_c(z) v_w = −∏_i(w₁−z_i)/∏_{m≤3}(w₄−w_m) (= 232551/1558000 at the point). 2 sites exact:
  ⟨0_w|t₁ᶻ(θ₂) = hθ₁₂(v_w⊗u_wG) − h²q⟨S|, ⟨0_w|t₁ᶻ(θ₁) = hθ₁₂(v_wGᵀ⊗u_w) − h²q⟨S|; inverses always localise
  (⟨0_w|K = 0 ⇒ resolvent trivial): t₁(θ₁)⁻¹ → vG^{−T}, t₁(θ₂)⁻¹ → uG⁻¹, t₃(θ₁−h) → vG^{−T}, t₃(θ₂+3h) → uG⁻¹,
  t₁(θ₁)t₁(θ₂)⁻¹ = Gᵀ⊗G⁻¹.
- Generic w (q ≠ 0): only end sites + end-growing chains localise (fixed G^{∓T}); interior single words entangle the end
  block; w₄ ∈ {z} changes nothing. Locus w₁ ∈ {z}: A sites are G_z-eigenvectors (frozen, all words scalars), F sites free:
  t₁(θ₂) → uG, t₁(θ₄)⁻¹ → uG⁻¹, t₁(θ₂)t₁(θ₄) → both by G; opposite powers entangle sites 2,4. One-parameter family
  ⟨0_w|(1 + a t₁(θ₂)) = ⟨0_w|ρ₂(1−κG), a = −κ/λ, λ = hθ₁₂∏_{k≥1}θ_{2,2k+1}(θ_{2,2k+2}−h) [exact L = 2,3 pairs]. Two-pair
  common family: ρ₂(g)ρ₄(g), g = 1−κG, in span{1,t₁(θ₂),t₁(θ₄),t₁(θ₂)t₁(θ₄)}: a = −κ(θ₂₄−h)/(θ₂₄λ), b = −κ/α₄,
  c = κ²/γ with ⟨0|t₁(θ₄) = α₄(ρ₄G + (h/θ₂₄)ρ₂G), α₄ = −hθ₄₁θ₄₂θ₄₃θ₄₅(θ₄₆−h), ⟨0|t₁(θ₂)t₁(θ₄) = γρ₂Gρ₄G,
  γ = λα₄(θ₂₄+h)/θ₂₄; three pairs: 8 words, d_S κ^{|S|}. Pairings v₃·u₂ = v₂·u₃ = 0 on the locus (extremal).
- Generic w word span on pair 1 (`twisted-vacuum-wordspan-L2.wb`, ≤ 2 factors of t₁,t₃): dim 𝒮 = 37, dim 𝒮₁₂ = 4 =
  span{v⊗u, vG^{−T}⊗u, vG^{−T}⊗uG, S − κv⊗uG}; rank-1 lines v(α+βG^{−T})⊗u and vG^{−T}⊗u(α+βG) only.
- Diagonal frame + untwisting (`twisted-vacuum-untwisting.wb`): locus → pair e₁⊗(e₄+αe₂+βe₃) (A frozen), overlaps only
  the BPS vacuum (sector "all A = Y¹"); generic w with w₂,w₃ = z₂,z₃ → e₁⊗(e₄+α̃e₁), sees [L−1,0,L−1] only. Neither is
  the note's W₋₁ (needs an O(1) A-site rotation).

## 2. 2603 reading (`notes/rv-sov-dictionary.md` §5; e-print in scratchpad)

Site = 6 = ∧²C⁴, vacuum Z = e₁∧e₂; κ-family = one-point word 1 − κt₊(θ₁) + κ²t₀(θ₁) (t₀ = ∧² transfer matrix) = group
element ∧²(1−κg); localisation (App. K) needs only ⟨Z|L_jk = 0 (j<k) + bonus L₂₁ = L₄₃ = 0, holds for any global rotation,
NO goodness condition; pairing of rotated with unrotated is quadratic in D (finite, non-extremal). Only one κ per block.
Untwisting: z, κ, ω → 1 together, rows of SoV matrices vanish, compensated by 𝒩_ℓ. Final C = 𝒩·𝕎_ω·𝔸_{ℓ,κ}/√𝔹.

## 3. Non-good reference W_a^L: exact families and continuous limits (`abjm-nongood-reference-limit.wb`, `odd-bridge-general-L.wb`)

- On ⟨Ω| = (v⊗u)^L with v·u = 0 (same pair everywhere; fails for independent pairs): P_{k−1}t₁(θ_k)^{∓1} = d_kP_k
  (A: inverse letters → Λ⁻¹, F: forward → Λ), 1/d_k = h∏_{s≠k}ϑ_ks (k odd; ϑ_ks = θ_ks for F, θ_ks+h for A),
  d_k = −h∏_{s≠k}ϑ̃_ks (k even; ϑ̃_ks = θ_ks for A, θ_ks−h for F); c_k = d₁⋯d_k [proved; exact L = 2,3; 1e-16 L = 4].
  c₁ = §5.5's formula; c₁₂ = c₂ = ∏_{A,s≥3}θ_{2s}/(θ_{1s}+h)·∏_{F,s≥4}(θ_{2s}−h)/θ_{1s};
  c_{2L} = ∏_{F<F'}(h²−θ²_{FF'})/∏_{A<A'}(h²−θ²_{AA'}) (L = 2: c₄ = −1108723/174000, c₃ = −321399232/112262625).
- Family (b) ⟨Ω|(α+βt₁(θ₁)⁻¹) = ⟨Ω|ρ₁(α+βc₁Λ⁻¹); family (a) ⟨Ω|t₁(θ₁)⁻¹(α+βt₁(θ₂)) = ⟨Ω|ρ₁(Λ⁻¹)ρ₂(αc₁+βc₁₂Λ) [exact].
  Both at every placement k = 1..2L via the chain. Products of (a),(b) are not product states.
- Untwisting (z = 1+εφ, θ = θ_ABJM+εθ⁰, ε = 10⁻²..10⁻⁴): (a) with α′+β′ = cε → bra (ε(c+φ₄)/z₁)N_ε,
  N_ε = ((e₁+a(z₁/z₄)e₄)⊗(e₄−a′e₁))⊗W_a^{L−1}, a′ = a(c+φ₁)/(c+φ₄); finite shift e^{−aE₁₄} → vac^L with one F site = e₁;
  ⟨N_ε|Ψ_ε⟩ → V₀ + (a−a′)⟨e₁e₁⊗vac^{L−1}|Ψ₀⟩ at rate O(ε), φ- and c-independent; ⟨e₁e₁…⟩ = ⟨e₄e₄…⟩ = C₁ (data):
  L = 2 [1,0,1]: 0.144338 ∓ 0.25i; L = 3 [2,0,2]: ±i/√24, 1/√24. (b) mirror with a″ = a(c+φ₄)/(c+φ₁); κ = z₄, a = −1
  reproduces §5.5 exactly at finite ε. Homogeneous ABJM-point limit REGULAR (c₁ = (2h²)^{−L}, c₁₂ = 1); singular locus
  θ_A = θ_F; normalised τ̂_k = d_k^{±1}τ(θ_k) finite (twisted vacuum: τ̂_odd = z₁, τ̂_even = z₄);
  ⟨bra_ε|Ψ_ε⟩ = [∏_{j<k}τ̂_j^{(−1)^j}](cε−1+τ̂_k)⟨W_a^L|Ψ_ε⟩. Even bridge not reached by these.

## 4. THE PAIR κ-FAMILY (`pair-group-rotation.wb`, `checks/session-2026-10-12/pairgroup/`)

- Pair-vacuum pattern (sympy): ⟨e₁e₄|(L_AL_F)_{jk}(u) = (u−θ_A)(u−θ_F)δ_jk⟨e₁e₄| − h(u−θ_A)δ_{k4}⟨e₁e_j| + h(u−θ_F)δ_{j1}⟨e_ke₄|
  − h²δ_{j1}δ_{k4}⟨S|: 7 identical zeros (lower + (2,3)), bonus zeros at u = θ_A ((j,4), j = 2,3) and θ_F ((1,k), k = 2,3).
- Lemma (2603 App. K style): rows with 1 ∉ I, 4 ∈ I are diagonal on vacuum sites; t₁(θ_F) and t̂₃(θ_A−h) localise on pair 1
  for any global rotation and twist; t₂(θ_F+h) localises ONLY at θ_A − θ_F = 2h (scan d = −5..5). Localised algebra on
  pair 1 at the ABJM point is 5-dim {1, X₁ = t₁(θ₂), X₂ = t₂(θ₂+h), X₃ = t̂₃(θ₂+h), X₃X₁} (4-dim off the point).
- On ⟨W_a^L| (pair-periodic θ, θ_A − θ_F = 2h, twist z₂ = z₃ =: y, s = 1−κ(z₁+z₄), λ_k(u) = ∏_{n≥2}∏_{j<k}(u−θ_{A,n}−jh)
  ·(u−θ_{F,n}−kh)∏_{j≤k−2}(u−θ_{F,n}−jh)):
  K_κ = s² + κs[X₁/(h²λ₁) + X₂/(2h⁴yλ₂) + X₃/(6h⁵y²λ₃)] + κ²X₃X₁/(24h⁷y²λ₁λ₃),
  ⟨W_a^L|K_κ = (1−κz₁)²⟨W_{a′}⊗W_a^{L−1}|, a′ = a(1−κz₄)/(1−κz₁)  [exact L = 2, 3; a ∈ {1/3,−1}; κ ∈ {3/10,−7/5,1/z₄}].
  Unique, a-independent. z₂ ≠ z₃: leftovers ∝ E₂₂+E₃₃ and z₃E₂₂+z₂E₃₃ (6 directions vs 5 operators) — fails; fails on
  generic BPS references (entangled, 3.5% off). Unstructured sums reach the target at L = 2 (meaningless) but not at
  L = 3 for z₂ ≠ z₃.
- One-pair even bridge κ = 1/z₄: ⟨vac⊗W_a^{L−1}|Ψ⟩ = (z₁²−z₁Σ+Π)/(z₁−z₄)²⟨W_a^L|Ψ⟩ (Σ, Π eigenvalue combinations, in nb).
- Untwisting (z = 1+ε(0.9,.35,.35,−1.1), pair-preserving inhomogeneity ε(0,.37,−.52)): (i) κ = 1/z₄, a = −1 → ⟨vac⊗W₋₁^{L−1}|Ψ₀⟩,
  O(ε), L = 2 singlets −1/√40, 1/√24 (= ABJM-md §5.7); (ii) κ → 1 with a′−a = −1 → ⟨vac^{L−1}⊗W₋₁|Ψ₀⟩ at L = 3, adjoints
  −0.153484 (×2), −0.035872 ∓ 0.028637i, −0.171126, −0.027138 (max err 3.8e-5 at ε = 1e-4). [num]

## 5. SoV dictionary at 4 sites (`sov-Q-dictionary-L4.wb`, `checks/session-2026-10-12/sovq/sov_l4.py`)

- Word basis ⟨a₁,k₂,a₃,k₄| = ⟨0|t_{a₁}(θ₁−h)t₁(θ₂)^{k₂}t_{a₃}(θ₃−h)t₁(θ₄)^{k₄}: rank 256, exactly z-independent, block-diagonal
  in the charge |x| = a₁+k₂+a₃+k₄ (= Σ(i−1)_A + Σ(4−j)_F; blocks 1,4,10,20,31,40,44,40,31,20,10,4,1). RV normalisation
  N_x = 1 in the Γ[V₄] gauge (polynomial gauge: N_x = ∏_F V₄(θ_F)^k∏_A∏_{b≤a}V₄(θ_A−bh); V₄(θ₂) = −2263/14000).
- Q convention Σ_{a=0}^{4}(−1)^aτ_a(u)Q(u−ah) = 0 (τ₄ = det G·qdet; NOTE shift Q^{RV}(u) = Q(u−h)); Q_i = z_i^{u/h}Φ(u)q_i,
  Φ(u) = V₄(u)Φ(u−h), q_i monic, deg q_i = 2+μ_i. Dictionary [hp 1e-42]: F sites τ₁(θ_F) = Q_j(θ_F)/Q_j(θ_F−h) for ANY j;
  A sites τ_a(θ_A−h) = D_a/D₀, D_a = det_{i∈I, j∈{1,2,3,4}∖{a+1}}Q_i(θ_A−jh) for ANY triple I. Separated variables (B):
  F: x ∈ {θ_F, θ_F+h} (k at θ_F+h); A: {θ_A,…,θ_A−3h}∖{θ_A−ah}. Full RV2 (5.20) product agrees.
- Dual ket basis: |0̄⟩ = e₄⊗e₁⊗e₄⊗e₁, A sites t₁(θ_A)^m, F sites t_n(θ_F+nh); ket wave functions mirror (A: powers of one
  dual Q; F: 3×3 det of three dual Q's). Gram ⟨x|y⟩ nonzero only if S_x+S_y = 12, triangular, 1456 nonzeros; NO diagonal
  measure; M = G⁻¹; completeness verified; 1/𝒩 = ⟨Ψᴸ|Ψ⟩/(⟨0|Ψ⟩⟨Ψᴸ|0̄⟩) = Σ τ̃_yM_{yx}τ_x. ⟨0|Ψ⟩ ≠ 0 for all 256 states.
  GLMR-type 12×12 determinant for 1/𝒩 does not fit. Calibration (i) holds only for global BPS rotations v⊗w⊗v⊗w.

## 6. One pair exact (`sov-onepair-determinant.wb`, `checks/session-2026-10-12/onepair/`)

- Basis graded (e_i⊗e_j grade i−j+3; blocks 1,2,3,4,3,2,1); c_{ak}(v,u) bilinear on the anti-diagonal i−j = a+k−3, e.g.
  c₀₀ = v₁u₄, c₀₁ = [(h−θ)v₁u₃ + hv₂u₄]/(hθ(h−θ)), c₀₃ = −v₁u₁/(h³θ²(h−θ)), c₃₃ = −v₄u₁/(6h⁶θ³(h−θ)(2h−θ)(4h−θ)) (θ = θ₁₂).
  Sparsity: generic BPS 16; one site rotated 7; A twist-adapted deg 1: 3; F: 2; both: 5.
- Determinant [exact 16/16]: ⟨v⊗u|Ψ⟩ = ⟨0|Ψ⟩det[w; Q₁(x_b); Q₂(x_b); Q₃(x_b)]/Q₁₂₃(θ₁−h), x_b = θ₁−bh, w_b = (−1)^bΣ_k c_{bk}r^k,
  r = τ₁(θ₂); closed form of the weight row: Num_b = u(R)[(bh−θ)v_{b+1} + hΣ_{i>b}v_iR^{i−b−1}] − hR^{3−b}𝒮(R), R = r/(hθ),
  u(R) = Σu_jR^{4−j}, 𝒮 = Σ_{i≥j}v_iu_jR^{i−j}; on product eigenstates overlap = u(R) × (A-site 4×4 det); zero-weight states
  add h𝒯(R)Cas(R^{u/h},Q̃₁,Q̃₂,Q̃₃). No single Q-only determinant (weight row cubic in r). Torus: ρ(X), X ∈ C[Gᵀ], equals
  Σc_{ak}(e₁X,e₄X^{−T})t_a(θ₁−h)t₁(θ₂)^k as an operator.
- Two pairs: 16×16 coefficient matrix rank 16 generically (3 if pair 2 is the vacuum; 2 for F-only deg 1); no normalisation
  factorises; ⟨0|t₁(θ₄) = hθ₁₄θ₃₄[h e₁₃₁₄ + θ₂₄e₁₄₁₃]; F-only: ⟨e₁u₂e₁u₄|Ψ⟩/⟨0|Ψ⟩ = (γ₂+δ₂ρ₂)(γ₄+δ₄ρ₄) − hγ₂δ₄ρ₂/θ₂₄
  − hδ₂δ₄ρ₂ρ₄/(θ₂₄+h), ρ₂ = τ₁(θ₂)/(hθ₁₂θ₂₃(θ₂₄−h)), ρ₄ = τ₁(θ₄)/(hθ₁₄θ₂₄θ₃₄). NOT a Vandermonde measure (the γ₂δ₄ cross
  term carries ρ₂, not ρ₄): det of the 2×2 bilinear form = hρ₂[ρ₂/θ₂₄ − ρ₄/(θ₂₄+h)] is not separable, so no one-row-per-site
  determinant exists in the left-basis/τ language (same for A–A: det = (h/θ₁₃)(1−ô₁ô₃)); A–F pairs have NO cross factor
  (det = 0) in the simple directions, but other directions do mix (site 1 + e₂ involves τ₁(θ₂), τ₁(θ₄), τ₁(θ₁−h)).
  A 2603-type determinant needs a C-type right SoV basis / non-diagonal measure for the mixed chain (open).

## 7. Two-block bras in the SoV basis (`sov-twoblock-expansion-L4.wb`, `checks/session-2026-10-12/twoblock/`)

- Operator-Schmidt rank of M⁻¹: 169 across {1,2}|{3,4}. Left block cheap: generic block {1,2} with {3,4} at the vacuum
  → 36 words (a₃ = 0, k₄ ≤ 2), rank 3; the reverse → 130 words, rank 16; common rotation 1−κG_cᵀ on pair 1 → 20 words.
- ⟨0|_comp = λ₁²(vac⊗vac)ρ(n′D′), λ_j = 1/∏_{m≠j}(z_j−z_m), r = λ/λ₁, n′ = ((1,1,1,r₄),(0,1,0,−r₂),(0,0,1,−r₃),(0,0,0,1)),
  D′ = diag(1,r₂,r₃,1) — finite as z → 1: the companion SoV vacuum is in the HW orbit of the physical vacuum. Block
  rotations giving the physical bra: g₁ = 1+b𝔼, g₂ = 1+b′𝔼 (𝔼 = E₁₄ in the twist eigenbasis); sparsest b′ = 0:
  B(b) = (e₁+bλ₁v_{z₄})⊗(e₄−bℓ₁)⊗e₁⊗e₄ (v_{z₄} = e₁∏_{m≠4}(G_cᵀ−z_m), ℓ₁ = λ₁e₄∏_{m≠1}(G_c−z_m)), 36 words (1/17/36 at
  b⁰/b¹/b²). Untwisting: ⟨B(b)|Ψ_ε⟩/(λ₁²(−b)^{2−w₁}) → data C₂ (singlets −0.15807 → −1/√40, 0.20408 → 1/√24) linearly in ε.
  36-term SoV sum = direct overlap to 1e-44 (exact c_x + 60-digit); double precision fails (cond M ~ 1e18–1e25).
- K_b = ⟨B(b)|Ψ⟩/⟨0|Ψ⟩ = Σ_aτ_a(θ₁−h)P_a(r₂,r₄) = det[Q_i(θ₁−βh); (−1)^{β+1}P_β(r₂,r₄)]/Q₁₂₃(θ₁−h): one bordered 4×4
  Casoratian; no finer structure (no 2×2 within a pair, F polynomials rank 3, CP rank 7, no sparser word basis).

## 8. Polarisation forms (`twoblock-polarisation-forms.wb`, `data/abjm-twoblock-forms.json`, `checks/twoblock-polforms/`)

- F_Ψ(v₂,u₂;v₃,u₃;M) = ⟨(v₂⊗u₂)^{L−M}⊗(v₃⊗u₃)^M|Ψ₀⟩ is a polynomial in the 8 Borel invariants a₂ = (v₂)₁, a₃, d₂ = (u₂)₄,
  d₃, mv = (v₂)₁(v₃)₂−(v₂)₂(v₃)₁, mu = (u₂)₃(u₃)₄−(u₂)₄(u₃)₃, p₂₃ = v₂·u₃, p₃₂ = v₃·u₂ (plus p₂₂,p₃₃ off-BPS) for all 184
  states; minimal and complete (7248 sectors). Exponents: mv^{e₅}mu^{e₆}, e₅ = (2q+r−p)/4, e₆ = (p+2q−r)/4;
  e₈+e₉ = s = L−(p+2q+r)/2; number of structures N = max(0, min(s,n₂−e₅,n₃−e₆) − max(0,s−n₃+e₅,s−n₂+e₆)+1) (0 mismatches/882).
  F/F_vac = Σc x^{e₈}y^{e₉}μv^{e₅}μu^{e₆}, x = p₂₃/(a₂d₃), y = p₃₂/(a₃d₂); note's C_M = F/F_vac at (1,−1,0,0), reproduced to 7e-15.
  Cyclicity c(M;e) = (−1)^{e₅+e₆}c(L−M; swapped); p = r states symmetric x ↔ y; one constant per (state, M) except
  [2,0,2] L = 4 M = 2 (two: C = 2c₀−c₁).
- The note's polarisation (μv = μu = 0 ⇒ x = −y under BPS) is BLIND to states with (x+y) | F/F_vac ⇔ s odd ⇔ L+p odd (the
  even-bridge parity rule is a property of the polarisation) and to all q > 0 states (need μv or μu ≠ 0); example P*:
  v₃ = (1,1,0,1), u₃ = (1,1,0,−2) gives F/F_vac = (x+y)/√12 ≠ 0 on L = 2 [1,0,1]. [2,0,2] E = 16 at L = 3: F = 0 for all
  two-block polarisations. 5 accidental zeros (E = 16, 32).

## 9. ABJM-point limit of the SoV basis (`sov-abjm-point-limit.wb`, `checks/abjmlimit/`)

- Pair-periodic θ: generic words have rank 70 for every d (coincident inhomogeneities; t_a(θ₁−h) = t_a(θ₃−h) etc.). Pair-by-
  pair θ = (θ_F1+d, θ_F1, θ_F2+d, θ_F2): rank 256 generic; d = 2h: 81 = 9²; one pair det M(d) ∝ d²³(d−h)¹¹(d−2h)⁷(d−3h)³(d−4h)⁴.
- Exact identities at u₀ = θ_A−h = θ_F+h: t₃(u₀) = 0; t₂(u₀) = t₁(u₀)t₁(θ_F). On states τ₃(u₀) = 0, τ₂(u₀) = τ₁(u₀)τ_F
  (the B collisions; 13 eigenvalues per pair). Q-reading: A-site Casoratian columns at θ_A−h, θ_A−2h proportional
  (D₀ = D₃ = 0, D₂ = τ_FD₁); LR map t̄₁(u) ↔ t̄₃(u+h) sends t₁(θ_F) to t̄₃(θ_A−h) exactly at d = 2h.
- Regularised basis (= d → 2h limit of RV basis, Jordan basis of B: 10 simple + 3 J₂ per pair): per resonant A site
  O(a,k) = {1, t₁(u₀), W = t₂′(u₀) − t₁′(u₀)t₁(θ_F) (k ≤ 2; t₂(u₀) at k = 3), Z = t₃′(u₀) = g₃′t̂₃(u₀)}; words
  ⟨0|O₁t₁(θ₂)^{k₂}O₃t₁(θ₄)^{k₄}, rank 256; ⟨x|Ψ⟩ = ⟨0|Ψ⟩∏_pairs λ(a,k)τ_F^k, λ ∈ {1, x = τ₁(u₀), w = τ₂′(u₀)−τ₁′(u₀)τ_F
  (k ≤ 2) or xτ_F (k = 3), y = τ₃′(u₀)} [hp 3e-48]. x, w, y from the u-derivative of the twisted Baxter equation (3×3 Cramer
  from q_j, q_j′ at θ_F+{2h,h,0,−h,−2h}). Fixed-bra coefficients diverge like δ^{−#resonant pairs}, overlaps converge.
  Physical fully-periodic point: needs shift-derivatives up to order 7 (homogeneous limit open).

## 9b. Locality brute force at 3 pairs, generic θ and z (`pair-rotation-locality-L3.wb`, `checks/session-2026-10-12/locality/`)

Diagonal frame, reference W_a^L, words on pair 1 (48 nonvanishing factors t_a(θ_s+mh), s ∈ {1,2}, m = −4..4; products ≤ 4
factors; exact mod p + Fractions). Localised subspace dim 6 at L = 2 and 3. Targets T1 = W_{a′}⊗W_a^{L−1}, T2 = common
rotation: in the span at ≤ 2 factors for L = 2 but only at ≥ 4 factors for L = 3 (rank 563 of 588), dense (no sparse
solution: OMP, L1, exhaustive ≤ 5 words all fail at L = 2). Obstruction: T1, T2 ∈ span{ψ₀, ψ₀t₃(θ₁−h), ψ₀t₃(θ₁−h)t₁(θ₂),
ψ₀t₁(θ₂)} ⊕ C·(S⊗W_a^{L−1}) with the isolated pair singlet S⊗w₀ not in the word span (coefficients identical at L = 2, 3:
ψ₀-coefficient 9(8a′−7)²/169, singlet 35(3a′−1)(8a′−7)/806). Needed word length grows ~2(L−1). Local universal
identities confirmed: c₁; c₃ (t₃(θ₁−h) → ρ₁(Λ⁻¹), c₃ = −(3/4)detΛ∏_F(θ_{1s}−h)(θ_{1s}−2h)(θ_{1s}−4h)∏_{A,s≥3}(θ_{1s}−h)(θ_{1s}−2h)
(θ_{1s}−3h)); c₁₂; ψ₀t₁(θ₂) = α[ψ₀ρ₂(Λ) + (σ/α)S⊗w₀], σ/α = ha(z₁−z₄)/θ₁₂ (L-independent). Calibration (i) holds only for
ψ₀ = W_a^L. ⇒ consistent with §4: the few-word solution exists only at the ABJM point with z₂ = z₃.

## 9c. ℓ-pair generalisation of K_κ and all even bridges at L = 3 (`even-bridge-lpairs.wb`, `checks/session-2026-10-12/lpairs/`)

Point: ABJM with distinct θ_F = (−3/5, −11/8, −7/9, 5/6), z = (3/2, 5/4, 5/4, 4/7), a = 1/3, r = z₁/z₄; eng.py matches pgbase.
- K^(1) works at L = 2, 3, 4 [exact]; t₂(u₀) = t₁(u₀)t₁(θ_F), t₃(u₀) = 0 confirmed as operator identities.
- NO common-κ ℓ-pair word: K^(p) (X's at pair p) changes pairs 1..p (block localisation); the 5^ℓ multilinear monomials
  block-localise but the target W_{a′}^ℓ⊗W_a^{L−ℓ} is not in their span (ranks 25 → 26 at (L,ℓ) = (2,2),(3,2),(4,2); 125 → 126
  at (3,3)); degree 3 spans the whole 285-dim cyclic module (trivial). The structure is PRETWISTED:
  ⟨W_{ar}^{p−1}⊗W_a^{L−p+1}|K^(p)_κ = (1−κz₁)²⟨W_{ar}^{p−1}⊗W_{a′}⊗W_a^{L−p}| (L = 3, 4, all p); (X₃X₁)^(p) = K^(p) at
  κ₀ = 1/(z₁+z₄) (s = 0); ∏_pK^(p)_{κ_p} works only with κ₁ = … = κ_{ℓ−1} = κ₀ — the end-growing chain, not a common block.
  Right-end words Y₁ = t₁(θ_{A,L}), Y₂ = t₂(θ_{A,L}), Y₃ = t̂₃(θ_{A,L}+h), Y₁Y₃ rotate pair L freely; left × right do not combine.
  ℓ = L = global ρ(1−κΛ), eigenvalue [(1−κz₄)/(1−κz₁)]^L∏(1−κz_i)^{μ_i}.
- Even bridges at L = 3 [f64, all 26 HW states, O(ε)]: finite-twist identity ⟨W_{a′}⊗W_a^{L−1}|Ψ⟩ = K_κ(τ)⟨W_a^L|Ψ⟩/(1−κz₁)²
  exact (1.3e-6 cancellation). Families: F1 (a′−a = −1, a = 1/3) and F4 (κ = 1/z₄, a = 1) → ⟨W₋₁⊗vac²| = C₄ (M = 1);
  F2 (a′−a = +1) and F3 (κ = 1/z₄, a = −1) → ⟨vac⊗W₋₁²| = C₂ (M = 2); G global → ⟨W₋₁³| = C₀. Symmetry C_M = (−1)^{L+μ₄}C_{L−M}
  for HW states ⇒ the one-pair word covers M = 1, L−1, L at every L; **first open case L = 4, M = 2.**
- Linear prescription κ = 1−cε with a = (φ₁−c)/(φ₁−φ₄); direction- and c-independent limits. z₂ = z₃ is a genuine
  finite-twist restriction (residual linear in z₂−z₃; symmetric choice y → (z₂+z₃)/2, y² → z₂z₃ is the LSQ optimum) but
  only z₂−z₃ = O(ε) is needed in the limit (identity breaks by ~0.1(z₂−z₃)²/ε; limits unchanged).
- Q-form (sov_l4 conventions, u₀ = θ_F+h, q_I = P_I/F₃ reduced Casoratian): τ_F = Q_j(θ_F)/Q_j(θ_F−h) (any j);
  x = τ₁(u₀) = [V₄^{rest}(u₀)/F₃^{rest}(θ_F)]·Z_I⁴·det[z⁻¹q(θ_F+h), z⁻³q(θ_F−h), z⁻⁴q(θ_F−2h)]/q_I(θ_F); τ₂(u₀) = xτ_F;
  τ̂₃(u₀) = h V₄^{rest}(u₀)V₄(θ_F)V₄(θ_F−h) Z_I [F₃^{rest}(u₀)/F₃^{rest}(θ_F)] q_I(θ_F+h)/q_I(θ_F) (any triple I) — a single
  dual-Q ratio. Checked L = 2 (256 states × 4 triples, 3.5e-8), L = 3 (33 states, 1e-6). K(κ) does not factor into
  rational linear factors in κ. Open: ⟨W_a^L|Ψ⟩ in SoV/Q form; middle M; norms.

## 9d. K_κ on the companion SoV vacuum; SoV-normalised formula (`even-bridge-sov-formula.wb`, `checks/session-2026-10-12/sovformula/`)

- On ⟨0|_comp the pair-1 localised span is 5-dim and z-independent: 1 → e₁e₄, X₁ → 2h²e₁e₃, X₂ = t₁(u₀)X₁ → 2h⁴(e₁e₂+e₂e₃),
  Z = X₃ → −12h⁵e₄e₄, X₃X₁ → −24h⁷e₄e₃; each is a single regularised SoV word (⟨0|, ⟨0|t₁(θ_F), ⟨0|t₁(u₀)t₁(θ_F), ⟨0|Z,
  ⟨0|Zt₁(θ_F)). ⟨0|K_κ has Schmidt rank 3 (not a product). Product states in the span: (ae₁+be₄)⊗(ce₃+de₄), BPS iff bd = 0.
- ⟨0|_diag = λ₁^L⟨vac^L|ρ(n′D′) exactly at every L (r_j ε-independent). No localised word puts e₁ on the F site ⇒ p₂₃ ≡ 0:
  the SINGLETS (form c·xy) are invisible from ⟨0| with pair-1 words; [1,0,1] only via y (odd bridge). A common rotation
  1−κΛ would not help: (n′f(Λ)n′⁻¹)₁₄ = f[z₁..z₄]/λ₁ = 0 for deg f ≤ 2 (the first-order E₁₄ entry (φ₁−φ₄)r₄ cancels since
  Σ_j r_jφ_j = 0). What reaches the singlets: 2603-type twist-dressed pair vectors, A by (ζ−Ĝᵀ)⁻¹, F by (ζ′−Ĝ)⁻¹, ζ,ζ′ = 1+O(ε),
  BPS iff ∏(ζ−z_j) = ∏(ζ′−z_j), all Borel invariants O(1) — not in the localised span, would have to be the reference.
- SoV-normalised formula [L = 2, mpmath 80 digits]: for any word K, |C_K(Ψ)| = |P_K(τ_Ψ)/P_K(τ_vac)|·|λ^{w−w_vac}𝒩(Ψ)/𝒩(vac)|^{1/2},
  𝒩 = ⟨0|Ψ⟩⟨Ψᴸ|0̄⟩/⟨Ψᴸ|Ψ⟩ = 1/Στ̃_yM_{yx}τ_x, 𝒩(vac) = (λ₁λ₄)^L, λ^w = ∏λ_j^{w_j} (= 2603's C = 𝒩𝕎𝔸/√𝔹 with 𝕎 ≡ 1, 𝔹 = 1/𝒩);
  rests on |⟨0|Ψ⟩|² = |λ^w𝒩|‖Ψ‖² (ratio → 1 like ε for all 6 L = 2 states, 3 φ's). ⟨0|-based limits are NOT the naive
  HW values: A-family on [1,0,1] → (y₀/3)⟨e₂e₂e₁e₄|Ψ₀⟩ = ⅓e^{∓iπ/3}× naive (observed, not derived). L = 3 not done.
- Q-forms at the ABJM point: P_K = c₀ + c₁τ_F + c₂xτ_F + c₃y + c₄yτ_F; x, w, y from Q_i′(u₀) − xD_i + wQ_i(u₀−2h) − yQ_i(u₀−3h) = 0,
  D_i = Q_i′(u₀−h) − τ_FQ_i′(u₀−2h), any 3 of the 4 Q's (τ_F to 4.5e-49; x,w,y to 3e-27). Norm at the ABJM point: ket side
  mirror identities t₂(θ_A) = t₁(θ_A)t₁(u₀), t₃(θ_A+h) = 0; regularised ket basis {1, t₁(u₀), t₂′(θ_A) − t₁(θ_A)t₁′(u₀),
  t₃′(θ_A+h)} (rank 256); Gram with the regularised bra basis: 3323 nonzeros, no longer triangular / charge-graded;
  1/𝒩 = Στ̃Mτ to 1.4e-5 (f64).
- CAVEAT: at z₂ = z₃ the companion frame is singular and the diagonal-twist chain has a gl(2) symmetry — no SoV basis;
  the K_κ route must be combined with SoV at z₂ − z₃ = O(ε) (sufficient for the limit per §9c).

## 10. Open at session end (see NEXT.md 2k–2m)

(all agents finished) Middle even bridges 2 ≤ M ≤ L−2 (first: L = 4, M = 2) — need a second integrable family or a
twist-dressed reference (§9d); ⟨W_a^L|Ψ⟩ in SoV/Q form (2603's 𝕎) and the dressed-reference overlap as a Q-bracket
determinant; the right SoV basis / measure of the mixed 4/4̄ chain (no diagonal measure found); the homogeneous
(fully pair-periodic) limit; z₂ ≠ z₃ at finite twist; the ⅓e^{∓iπ/3} factor of §9d; L = 3 untwisting of the SoV formula;
rerun all notebooks after the licence renewal.
