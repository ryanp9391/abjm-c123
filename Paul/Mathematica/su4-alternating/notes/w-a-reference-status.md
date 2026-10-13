# What ⟨W_a^L| is, and how its overlaps become twisted-vacuum overlaps (2026-10-13)

For: Paul. Answers the confusion about ⟨W_a^L| in `plan-2026-10-13.md` §0: which object it is, in which frame the Ψ live, and
whether its overlap is in Q's. Conventions as in the README (A F A F, bras row vectors, A: v → v g, F: u → u g^{−T}).
Agent reports of the day: `agent-2026-10-13-reference-map.md` (A), `agent-2026-10-13-Wa-sov.md` (B),
`agent-2026-10-13-epsilon-expansion.md` (C).

## Headline

⟨W_a^L| is **not** itself a twisted vacuum: no twist sharing the companion SoV basis has it as an eigenbra (§3). But on every
state it can see, its overlap coincides with that of the companion twisted vacuum ⟨0_w| with **w = (w₁, z₂, z₃, w₄)**:

  ⟨0_w|Ψ⟩ = (v₁u₄)^L ⟨W_a^L|Ψ⟩  for every Ψ of weight (μ,0,0,−μ),
  a = v₄/v₁ = −(z₁−z₂)(z₁−z₃)(w₁−z₁) / ((z₄−z₂)(z₄−z₃)(w₁−z₄)),

and ⟨W_a^L|Ψ⟩ = 0 for every other weight. [exact to 40 digits, L = 2 pairs, all eigenstates of the sectors (2,0,0,−2),
(1,0,0,−1), (0,0,0,0); the argument is weight counting and holds at any L.] The overlap ⟨0_w|Ψ⟩ is a finite SoV sum in
Q-functions (agent B, verified to 1e-44 on all 256 states at L = 2). Hence **⟨W_a^L|Ψ⟩ is a Q-expression.** Simplest form (w₁ → ∞, §5): the SoV vacuum with every F site rotated by
(G_c − z₂)(G_c − z₃), i.e. a polynomial in the F-site τ₁(θ_F); the L = 2 even bridge assembled from it reproduces the data
singlets 1/√40, 1/√24 within 1e-11 at the ABJM point, Q-complete up to the SoV norm ⟨0|Ψ⟩/‖Ψ‖.

So Paul's remark of 2026-10-12 ("⟨W_a| is the vacuum of a different twist; I know the overlap in Q's") was right in substance.
The confusion was about *which* twist. It is not the triangular conjugate G_a = e^{aE₄₁}Λe^{−aE₄₁}. It is the companion twist
G_c(w) whose eigenvalues share z₂, z₃ with the twist of Ψ, while w₁ encodes a (a Möbius map) and w₄ only normalises.

## 1. Frames

- Diagonal (physical) frame, twist Λ = diag(z). Per pair W_a = (e₁ + a e₄) ⊗ (e₄ − a e₁); ⟨W_a^L| = ⟨vac^L| ρ(e^{aE₁₄}),
  vac = e₁ ⊗ e₄: a global Borel rotation of the physical vacuum.
- **The Ψ in every even-bridge formula are eigenstates of t^Λ(u) in the diagonal frame.** K_κ is a word in t^Λ transfer
  matrices; the identity ⟨W_a^L|K_κ = (1−κz₁)²⟨W_{a′}⊗W_a^{L−1}| and the untwisting z → 1 are diagonal-frame statements.
- The companion frame enters only for SoV: |Ψ_c⟩ = R|Ψ⟩, ⟨B_c| = ⟨B|R^{−1}, R = V^{−T} on A ⊗ V on F, t^{G_c(z)}R = R t^Λ. Every
  bra is transported with the R of the twist z of Ψ, including ⟨0_w| (whose own twist is w). Eigenvalues and Q-functions are
  frame-independent, so P_κ(τ) is the same in both frames.
- ⟨W_a^L|Ψ⟩ = (aⁿ/n!)⟨vac^L|ρ(E₁₄)ⁿ|Ψ⟩ with n = L − μ₁ fixed by the weight (μ₁,0,0,−μ₁) of Ψ.

## 2. The four reference bras

All four are BPS product states (v ⊗ u)^L with v·u = 0. "Good" = the README criterion u Λ vᵀ = 0. λ_j = 1/∏_{m≠j}(z_j − z_m);
v(x) := (x³, x², x, 1); ũ_j := row j of V^{−1} (left eigenvectors of G_c(z)).

| bra | diagonal frame (per pair) | companion frame (per pair) | what it is | u Λ vᵀ | overlap with Ψ |
|---|---|---|---|---|---|
| ⟨vac^L\| | e₁ ⊗ e₄ | v(z₁) ⊗ ũ₄ | Bethe vacuum of t^Λ | 0 | 0 for Ψ ≠ vac |
| ⟨0\|_comp (SoV / GT vacuum) | λ₁ (λ/λ₁) ⊗ (1,1,1,1) | e₁ ⊗ e₄ | RV GT vacuum; z-independent word basis on it | 0 (good) | ⟨0\|Ψ⟩ = SoV normalisation; does not admit the five-word K_κ |
| ⟨0_w\| (twisted vacuum) | v_j = λ_j ∏_{m≠j}(w₁−z_m), u_j = ∏_{m≤3}(z_j−w_m)/∏_{m≤3}(w₄−w_m) | v(w₁) ⊗ u_w, u_w = row 4 of V(w)^{−1} | vacuum eigenbra of t^{G_c(w)}; same SoV basis as Ψ | −∏_i(w₁−z_i)/∏_{m≤3}(w₄−w_m) | finite Q-sum 𝕎_w (§4 B); admits no K_κ-family at z₂ ≠ z₃ (§4 A) |
| ⟨W_a^L\| | (e₁ + a e₄) ⊗ (e₄ − a e₁) | (v(z₁) + a v(z₄)) ⊗ (ũ₄ − a ũ₁) | vacuum eigenbra of G_a = Λ + a(z₁−z₄)E₄₁ (and of Λ(w) + a(w₁−w₄)E₄₁ for every w); never of a companion twist | a(z₄ − z₁) | = (v₁u₄)^{−L}⟨0_w\|Ψ⟩ with w = (w₁,z₂,z₃,·) on all sectors it sees; admits the five-word K_κ |

Remarks.
- The README goodness quantity u Λ vᵀ equals the companion-frame obstruction q(w,z) of the twisted-vacuum notebooks
  (232551/1558000 at the standard w, z).
- K_κ-admissibility at z₂ = z₃ (agent A, proven with `Reduce`): a common pair rotation by the five pair words exists iff either
  v, u are BPS inside every spectral block of Λ ({1}, {2,3}, {4}), or v, u are supported on indices {1,4} (torus images of W_a).
  The first class has vanishing cross-block pairings, so W_a is the only admissible reference that can see the bridge.
  The five-word coefficients of K_κ are specific to the W class; the block-BPS class has its own class-dependent words.

## 3. Why ⟨W_a^L| is not a twisted vacuum, and why that does not matter

**The bra is not.** Its transfer-matrix eigenvalue is τ_vac (conjugating the twist does not change eigenvalues), so its
Q-functions are those of ⟨vac^L|. If a twist sharing the companion SoV basis had ⟨W_a^L| as eigenbra, its dual-basis wave
functions would equal those of ⟨vac^L|, and completeness would force ⟨W_a^L| ∝ ⟨vac^L|. Agent B confirms this directly: the
SoV coefficients of ⟨W_a^L| are dense (all 256 words, every power aⁿ), and ⟨W_a|y⟩ does not factorise.

**The overlaps are.** Take w₂ = z₂, w₃ = z₃ in the twisted vacuum. Then u₂ = u₃ = 0: the F sites carry only e₁ and e₄. An A site
in e₂ or e₃ would then give the bra a positive 2- or 3-weight that no F site cancels. So on weights (μ,0,0,−μ) only the A
components v₁, v₄ contribute, each pair acts as (v₁e₁ + v₄e₄) ⊗ (u₁e₁ + u₄e₄) = v₁u₄ W_a, and v·u = 0 gives a = v₄/v₁ = −u₁/u₄.
The A sites of ⟨0_w| still have e₂, e₃ components, which is why the two bras differ while their relevant overlaps agree.
The map a(w₁) sends w₁ = z₁ to a = 0 (the vacuum) and w₁ = z₄ to a = ∞ (the lowest-weight vacuum e₄ ⊗ e₁).

## 4. Results of the three routes (agents of 2026-10-13)

- **A. Changing the reference to ⟨0_w| — negative** (`reference-admissibility-Kkappa.wb`). At z₂ ≠ z₃ no K_κ-family exists on
  ⟨0_w| unless its A sites freeze (w₁ ∈ {z}); eleven directions are missing from the five-word span, so a sixth localised
  operator cannot help. At z₂ = z₃ the diagonal image of ⟨0_w| is dominated by the A site (0,1,−1,0) (λ₂, λ₃ diverge); it
  then admits an F-site-only rotation K^F = 1 − κe₁′ + κ²e₂′ + 2κX₁ + 8κ²X₂ [exact, L = 2, 3], which is not the W_a family.
  GL(4) covariance: ρ(g)t^Λρ(g)^{−1} = t^{G′}, G′ = (gΛg^{−1})ᵀ. Untwisting of a rotated pair uses D = diag(c − φ_j).
  This route is superseded by §3: we keep W_a for the κ-family and use ⟨0_w| only for the overlap.
- **B. SoV expansion** (`Wa-sov-expansion-L4.wb`). ⟨W_a^L|R^{−1} = [(v(z₁)+a v(z₄)) ⊗ (ũ₄−a ũ₁)]^{⊗2} exactly. The bra is dense in
  the SoV words, and no determinant form exists for Σ_x c_x τ_x: the A–A coefficient matrix has rank 4, and products or sums
  of bordered Casoratians fail. The twisted vacuum has factorised dual wave functions with vacuum dual-Q's of twist w, and
  **𝕎_w(Ψ) := ⟨0_w|Ψ⟩/(⟨0_w|0̄⟩⟨0|Ψ⟩) = Σ_{x,y} τ̃_y(vac_w) M_{yx} τ_x(Ψ)** matches the direct overlap on all 256 states to 1e-44.
  At w = z it vanishes on all 255 non-vacuum states (biorthogonality as a finite Q-sum).
- **C. ε-expansion** (`even-bridge-epsilon-expansion.wb`; the factorisation, step 0a, the T₁₄ identity and K₁ are ported [wb]: exact −1/√40, 1/√24 at L = 2, 1e-14 at L = 3; ε-scans remain python). ⟨W_a^L|Ψ_ε⟩ = (aⁿ/n!)εⁿω_n; the M = 1 even
  bridge sees only n = 2 states ([L−2,0,L−2]). It factorises:
  **C_{M=1}(Ψ₀) = K₁(τ_Ψ₀) · Ω(Ψ₀)**, with K₁ = P_{κ=1} at the untwisted ABJM point = 2 − A − B − C (A = τ_F/(h²λ₁),
  B = τ₂(u₀)/(2h⁴λ₂), C = τ̂₃(u₀)/(6h⁵λ₃), and AC = 4 on every state), and Ω = lim ⟨vac|ρ(E₁₄)²|Ψ_ε⟩/(2ε²(φ₁−φ₄)²) bridge
  independent [wb exact at L = 2, 1e-14 at L = 3; ε-limits num]. A·C = 4 because X̂₃X₁ ∝ the one-pair shift (exact at L = 2), as in U ∝ S of `abjm-bridges.wb`. Exact handle: t^{G_a} = t^Λ + a(z₁−z₄)T₁₄, so
  (τ_vac(u) − τ_Ψ(u))⟨W_a^L|Ψ⟩ = a(z₁−z₄)⟨W_a^L|T₁₄(u)|Ψ⟩. No even/odd ratio relation exists: odd one-field bridges see n = 1,
  the even M = 1 bridge sees n = 2.

## 5. The even bridge M = 1 in Q-functions — verified at L = 2 (agent B follow-up, `Wa-sov-expansion-L4.wb` cells 30–48)

**Simplest member of the family: w₁ → ∞.** Then a → κ₀ = −(z₁−z₂)(z₁−z₃)/((z₄−z₂)(z₄−z₃)) and the twisted vacuum becomes, in the
companion frame, (e₁ ⊗ ũ)^L with ũ = e₄ᵀ(G_c − z₂)(G_c − z₃) = (0, 1, −(z₂+z₃), z₂z₃). The A sites are exactly those of the SoV
vacuum ⟨0|; only the F sites are rotated, by the degree-2 twist polynomial f(x) = (x−z₂)(x−z₃). In the diagonal frame f(Λ)
kills e₂, e₃, which is the weight argument of §3 in one line. Hence, on the sectors (μ,0,0,−μ), n = L − μ:

  ⟨W_a^L|Ψ⟩ = (a/κ₀)ⁿ (λ₁u₄)^{−L} ⟨0| ∏_{F sites} ρ_F(f(G_c)) |Ψ⟩,   u₄ = (z₄−z₂)(z₄−z₃),

and the right-hand side is an F-only SoV sum: at L = 2, ⟨W_a|Ψ⟩/⟨0|Ψ⟩ = (a/κ₀)ⁿ ũ·P(ρ₂,ρ₄)·ũ/(λ₁u₄)², with P a 3×3 matrix of
bidegree (2,2) in the F-site single-Q ratios ρ_s ∝ τ₁(θ_s) = Q_j(θ_s)/Q_j(θ_s−h) (the F–F exchange terms make it non-separable;
entries in `agent-2026-10-13-Wa-sov.md`). [hp: 1e-51 on all 54 sector states, generic θ and pair-by-pair ABJM θ; coefficients
reproduce at a second (h,θ) point.] The general-w Q-sum 𝕎_w is a degree-8 polynomial in w₁ equal to this times an explicit
factor [1e-46]. The dense 256-word expansion of ⟨W_a| (route B) is the unique representation on all states because it also
encodes the off-sector zeros; on the sector states 9 words suffice.

**General L — factorised form (agent D, `Wa-Fonly-L6.wb`, `agent-2026-10-13-Fonly-L6.md`).** At L = 3 the bra ⟨(e₁⊗ũ)^{⊗3}| is
exactly a combination of the 27 F-only words ⟨0|∏_F t₁(θ_s)^{k_s}, k_s ≤ 2 [exact, all 4096 components]. The coefficients have a
closed form that FACTORISES over the two middle twist eigenvalues. With normalised F-site eigenvalues y_s := τ₁(θ_s)/κ_s
(κ_s explicit θ-products in the notes) and ν_st = θ_s − θ_t over F sites:

  ⟨W_a^L|Ψ⟩ / ⟨0|Ψ⟩ = (a/κ₀)ⁿ (λ₁u₄)^{−L} · B_{z₂}(y) · B_{z₃}(y),
  B_z(y) = Σ_{S ⊆ F sites} (−z)^{L−|S|} g_S ∏_{s∈S} y_s,   g_S = ∏_{s<t∈S} ν_st² / (ν_st² − h²).

B_z is multilinear, so the product has degree ≤ 2 per F site (the 3^L words); for h → 0 it is ∏_s(y_s − z). [exact on the full
chain at L = 2, 3 (three points incl. the ε-chain); model-exact at L = 4; general-L proof open.] Q-form on states at L = 3:
1e-51 on 9 eigenstates of all four sectors. In the untwisting limit z₂, z₃ → 1 the two factors coincide, 𝕎 → B₁(y)²: a perfect
square (agent E finds the F-site factor → 25, 9 on the two L = 2 singlets).

**Even bridge.** Combining with the K_κ identity:

  ⟨W_{a′} ⊗ W_a^{L−1}|Ψ⟩ = P_κ(τ_Ψ) · (a/κ₀)ⁿ (λ₁u₄)^{−L} ⟨0|Ψ⟩ · ũPũ(τ_F's) / (1−κz₁)².

At L = 2 on the ABJM ε-chain (pair-by-pair θ, z = 1+εφ, κ = 1−cε, a = (φ₁−c)/(φ₁−φ₄), 90-digit eigenvectors) this equals the
direct overlap to 1e-66 on all 28 zero-weight states, and the singlets extrapolate to 1/√40 and 1/√24 within 1e-11 (up to the
phase of Ψ). The regularised SoV bases of reports §9 are not needed: only the F-only words ⟨0|t₁(θ₂)^{k₂}t₁(θ₄)^{k₄}, k ≤ 2,
enter, and they stay regular at the ABJM point.

**The norm in Q's** (𝔹 is known to Paul; this check was not needed but confirms the assembly) **— L = 2 fully Q-based (agent E, `abjm-sov-norm-L4.wb`, `agent-2026-10-13-sov-norm.md`).** ⟨0|Ψ⟩/‖Ψ‖ is
replaced by 1/√𝔹 with 𝔹 = 1/𝒩 = Σ_{x,y} τ̃_y M_{yx} τ_x, built with the regularised bra and ket bases at the ABJM point from eigenvalue
data only: bra side τ_F, x = τ₁(u₀), w = τ₂′(u₀) − τ₁′(u₀)τ_F, y = τ₃′(u₀); ket side τ_A, w̃ = τ₂′(θ_A) − τ₁(θ_A)τ₁′(u₀),
ỹ = τ₃′(θ_A+h); M = Gram(θ)⁻¹ a θ-dependent number table. Final L = 2 formula:

  |C_{M=1}| = lim_{ε→0} 𝒩_pre · 𝔸 · 𝕎 / √𝔹,   𝔸 = P_κ,  𝕎 = ũPũ = B_{z₂}B_{z₃},  𝒩_pre = (a/κ₀)²/((λ₁u₄)²(1−κz₁)²) → 1.

Limits (𝔸, 𝕎, 𝔹) → (10, 25, 2 500 000) and (2, 9, 7776) on the two singlets, i.e. exactly 1/√40 and 1/√24; extrapolation over
ε = 1/50 … 1/1600 gives 4e-16 and 9e-17. [hp; independent python cross-check.] The formula gives |C|; the sign is the phase of Ψ.
Finite-twist subtlety: the exact Q-statement is ⟨0|Ψ⟩²/⟨Ψᴸ_*|Ψ⟩ = λ^w/Σ τ̃Mτ (left eigenvector normalised by ⟨Ψᴸ_*|0̄⟩ = ⟨0̄ᵀ|Ψ⟩);
the Euclidean ‖Ψ‖² differs by H = 1 + O(ε) because twisted eigenvectors are not orthogonal — the limit is unchanged. Observed, not
proved: on the 16 parity-symmetric zero-weight states the ket data equal the bra data of the other pair (τ_A1 = τ_F2, x₁ = x₂,
w̃₁ = −w₂, ỹ₁ = −y₂, …), so 𝔹 depends on 7 numbers.

**Analytic untwisted limit at L = 2 (agent F, `even-bridge-Qform.wb` Part 3, `agent-2026-10-13-Qform.md`).** With 𝔹 the known
SoV norm (Paul), the two bridge/overlap factors have closed limits in untwisted one-pair eigenvalues at the ABJM point, with
τ = τ₁(θ_F), τ′ = ∂_uτ₁(θ_F), x = τ₁(θ_F + h), h = ½:

  |C_{M=1}| = K₁ · B₁² / √𝔹₀,   K₁ = 2 − 8τ − 1/(2τ) + 64xτ  (translation-invariant states),   B₁ = 1 + 4τ + 2τ′ + 16τ².

- 𝔸 → K₁: needs the regularised t̂₃(u₀) = ½t₃″(u₀) (double zero of the t₃ prefactor at θ_F1 = θ_F2). Singlets: (A,B,C) = (2,−12,2) and
  (−2,4,−2), K₁ = 10, 2 [exact].
- 𝕎 → B₁²: the normalised y_s diverge like 1/ν (ν = θ_F1 − θ_F2), but y₂ + y₄ → −(τ + hτ′)/(4h⁴) and g·y₂y₄ → τ²/(16h⁸) are finite.
  Singlets τ = ±¼, τ′ = 1 give B₁ = 5, 3 [exact]; on all 28 zero-weight states the limit formula holds with O(ε) convergence [hp].
  τ′ is genuine extra Q-data (a parity shortcut B₁ = 16τ² − 4τ − 1 + 8x works only on the singlets).
- At L = 3 the factors are not separately finite (𝕎 ~ ε², ⟨0|Ψ⟩/‖Ψ‖ ~ ε⁻³, 𝒩_pre ~ ε): only the product is; the L ≥ 3 limit needs
  the combined scaling.

Structure: 2603's C = 𝒩·𝕎·𝔸/√𝔹 with 𝔸 = P_κ (pair word, Q's at one pair), 𝕎 = ũPũ (F-site τ₁(θ_F) polynomial with exchange
terms), √𝔹 = SoV norm. Not yet a determinant: the F–F exchange structure of P and the norm both need the mixed-chain measure.

L = 3 even bridge (agent D): the Q-form equals the direct value to 1e-41 at every ε, and the extrapolated |C| of the six
[1,0,1] states match the data to 8 digits (norm still from the eigenvector at L = 3; the L = 3 Gram matrix is 4096 × 4096, too heavy for the same hp approach).

Open: proof of the B_z form at general L (and a determinant/Pfaffian form for B_z); the norm at the ABJM point
in Q's; the middle bridges 2 ≤ M ≤ L−2; determinant forms; the homogeneous limit.
