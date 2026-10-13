# F-only SoV form of ⟨W_a^L|Ψ⟩ at L = 3 pairs (6 sites) — agent report, 2026-10-13

Notebook: `Wa-Fonly-L6.wb` (all computations in the shared wolfbook kernel, prefix `fd`, no memoisation, vectors only).
Python cross-checks of the *model* in §2 (symbolic/exact, sympy): `checks/session-2026-10-13/fonly/` (`model_cherednik.py`, `model_L3.py`,
`check_conj.py` (run `model_L3.py` first: it writes the pickle `check_conj.py` reads), `check_conj_L4.py`); python was used only for these model checks, written while the shared kernel was blocked.
Conventions exactly as in the task (chain A F A F A F, site 1 most significant, bras = row vectors acted on from the right, Lax with
ρ_A(E_ij) = e_ij, ρ_F(E_ij) = −e_ji, t₁^G = Σ T_ij G_ji, companion twist G_c(z), R = V^{−T} ⊗ V ⊗ …).
Tags: [exact] rational arithmetic; [hp] fixed-precision arbitrary precision (digits stated); [num] machine precision; [python-exact] sympy rational.

## Headline

1. **Yes: at L = 3 the bra ⟨(e₁⊗ũ)^{⊗3}| (companion frame) lies exactly in the span of the 27 F-only words
   ⟨0|t₁(θ₂)^{k₂}t₁(θ₄)^{k₄}t₁(θ₆)^{k₆}, k_s ≤ 2.** [exact] rank 27 = rank with the target, residual 0 on all 4096 components, generic point; the
   27 words are independent and supported on exactly the 27 components e₁⊗{e₂,e₃,e₄}^{⊗3}. Same pipeline at L = 2 reproduces the known 9-word identity.
   No extension (k_s = 3, A-site words) is needed.
2. **Closed form of the coefficients, any L (found in an F-chain model, verified exactly on the full chain at L = 3):**
   ⟨(e₁⊗ũ)^{⊗L}| = ⟨0| B_{z₂}(Ŷ) B_{z₃}(Ŷ) with the single-root global F-rotation operator
   **B_z(Ŷ) = Σ_{S⊆{F sites}} (−z)^{L−|S|} g_S ∏_{s∈S} Ŷ_s,  g_S = ∏_{s<t∈S} ν_st²/(ν_st² − h²),  ν_st = θ_s − θ_t,**
   and t₁(θ_s)⟨0| = κ_s ⟨0|Ŷ_s (κ_s below). So c_k = [Ŷ^k](B_{z₂}B_{z₃}) / ∏κ_s^{k_s}. It is *not* a product ∏_s(α_s+β_s t+γ_s t²): the
   pair exchange factors g_st − 1 = h²/(ν_st² − h²) enter at every level, and the c_k are not symmetric in the s (only covariant under
   simultaneous permutation of (θ_s, Ŷ_s)).
3. [hp] Q-form on states verified at the generic point, 9 eigenstates in the sectors n = 0..3: relative deviation ≤ 2.5e-51 (50-digit arithmetic).
4. Even bridge at L = 3: see §4 (table). C_Q equals C_direct to 1e-40 at every ε, and the ε → 0 extrapolation is compared with the data.

## 1. Bra identity at L = 3 [exact]

* Words W_k built by propagating the bra ⟨0| = (e₁⊗e₄)^{⊗3} site by site through t₁^{G_c(z)}(θ_{2s}) (bra action = transposed site operators); each
  application 0.3 s. Machinery tested against explicit 16×16 and 256×256 matrices (ket and bra) and the frame map t^{G_c}R = R t^Λ verified
  exactly on 1 and 2 pairs.
* A first run reported a spurious "28th component" and a negative result also at L = 2: `Position[v, Except[0], {1}]` returns the head
  position 0 unless `Heads -> False` is given (notebook cells 9–12 carry the fix; kept in a markdown note).
* Level grading: ⟨0|t^k lies in total level |k| = k₂+k₄+k₆ (level = number of twist raisings e₄→e₃→e₂ on the F sites), so the 27×27 system
  is block diagonal with block sizes 1,3,6,7,6,3,1.
* The exact coefficient vector at the generic point is `fdSol` in the notebook (27 rationals, ordering `fdWordsK` = Tuples[{0,1,2},3] for (k₂,k₄,k₆)).

## 2. Structure of the c_k

**F-chain model.** On bras whose F-site vectors have no e₁ component, the A–F exchange operators act as scalars, and
t₁(θ_s)⟨v| = κ_s ⟨v|Ŷ_s with the commuting (Cherednik-type) operators on levels (C³)^{⊗L} (e₄,e₃,e₂ = level 0,1,2, N = raising, e₂N = 0)

  Ŷ_s = R_{s,s+1}⋯R_{s,L} N_s R_{s,1}⋯R_{s,s−1},   R_{sf} = 1 − (h/(θ_s − θ_f)) P_{sf}   (row-vector action, product in this order).

κ_s = h ∏_{A sites a}(±(θ_a − θ_s)) ∏_{F sites f≠s}(±(θ_f − θ_s)) with + when the other site is to the left of s and the sign reversed
(θ_s − θ_a, θ_s − θ_f) when it is to the right. Checks:
* [exact] κ_s read off the exact L = 3 words equals this formula (θ-generic point); L = 2: κ₂ = hθ₁₂θ₂₃θ₂₄, κ₄ = hθ₁₄θ₂₄θ₃₄.
* [exact] model words κ^k⟨vac|Ŷ^k == exact chain words, all 27 words × 27 components, deviation 0 (notebook cell 29).
* The target is the tensor ⊗_s f(N_s) on the vacuum, f(x) = (x − z₂)(x − z₃) = p − s x + x² (p = z₂z₃, s = z₂+z₃).

**L = 2 (from the exact P of the L = 2 report).** With Ŷ₂ = ((ν−h)/ν) ρ₂, Ŷ₄ = ρ₄ (ν = θ₂₄; ρ's as in `agent-2026-10-13-Wa-sov.md`) the nine entries P_bb′
collapse to a symmetric polynomial: ũPũ = ∏_{z∈{z₂,z₃}} [ (Ŷ₂−z)(Ŷ₄−z) + (g−1)Ŷ₂Ŷ₄ ], g = ν²/(ν²−h²), i.e. B_z = z² − z(Ŷ₂+Ŷ₄) + gŶ₂Ŷ₄.
(Level by level: e₂-type terms get g, e₁ terms are unchanged; [python-exact] the model reproduces this identically.)

**L = 3 and general L.** Closed form of §Headline 2.
* [exact, full chain] `fdClosedSol === fdSol`: the formula reproduces all 27 exact coefficients at the generic point (notebook cell 31).
* [python-exact] model: L = 3 at two (h, θ, z) points (first with p, s symbolic), L = 4 at one point (81-dim): ⟨vac|B_{z₂}B_{z₃} equals the target tensor exactly.
* Level-2 example, coefficient of Ŷ_sŶ_t (third site at level 0): p[(s²−2p) g_st + 2p] — the pair factors are g_st = ν_st²/(ν_st²−h²) with different values
  for the three pairs (484/259, −16/425, 3364/2139 at θ_F = (1/3,−2/5,3/7), h = 1/2), so the c_k are not symmetric in k.
* Reading: B_z(Ŷ) is the operator ∏_s ρ_s(G − z) (global rotation of the F sites by one root) written in the Ŷ's on the cyclic module ⟨0|C[Ŷ]/(Ŷ_s³);
  the quadratic twist polynomial is just the product of the two root operators, which commute. g_S = ∏_{pairs∈S}(1 − h²/ν²)^{−1}: every pair inside S pays one factor,
  so the structure is a "Wick/Jastrow" cluster expansion with pair weight g_st only (no genuine three-body weights up to L = 4).
* Not proven: the formula for general L is established exactly only up to L = 4 in the model and L = 3 on the full chain; the proof (an identity in the
  degenerate affine Hecke algebra for e_k(N)-type clusters) and the full-chain check at L = 4 are open. Not tested: complex/degenerate θ configurations (e.g. coinciding θ_F, where g_st → 0/0 is controlled by ν → 0 of the
  ABJM limit, g_st → 0 as ν → 0 at fixed h).

## 3. Q-form on states, generic point [hp, 50 digits]

⟨W_a^3|Ψ⟩ (diagonal frame, a = 3/7) vs (a/κ₀)ⁿ(λ₁u₄)^{−3}⟨0|Ψ⟩ Σ_k c_k τ₁(θ₂)^{k₂}τ₁(θ₄)^{k₄}τ₁(θ₆)^{k₆}; Ψ = eigenvectors of t₁^Λ(1/3) in the weight sector
(3−n,0,0,n−3), machine eigensystem then 3 steps of fixed-precision inverse iteration (50 digits), τ's by vector application.

| n | sector | dim | # states | max relative deviation | |⟨W|Ψ⟩| range |
|---|---|---|---|---|---|
| 0 | (3,0,0,−3) | 1 | 1 | 1.5e-57 | 1 |
| 1 | (2,0,0,−2) | 24 | 3 | 1.8e-55 | 0.00276 – 0.0186 |
| 2 | (1,0,0,−1) | 141 | 3 | 1.9e-54 | 3.2e-4 – 4.4e-4 |
| 3 | (0,0,0,0) | 256 | 2 | 2.5e-51 | 3.8e-4 – 1.4e-2 |

(the deviations are at the working precision; no state fails). Sector dimensions 1, 24, 141, 256.

## 4. Even bridge M = 1 at L = 3, ABJM ε-chain

Set-up exactly as the L = 2 assembly (`Wa-sov-expansion-L4.wb` cells 43–45): z = 1 + εφ, φ = (9/10, 7/20, 3/5, −11/10) (φ₂ ≠ φ₃, y = (z₂+z₃)/2),
θ = (θ_F1+2h, θ_F1, θ_F2+2h, θ_F2, θ_F3+2h, θ_F3), θ_F = ε(0, 37/100, −52/100), c = 3/10, κ = 1 − cε, a = (φ₁−c)/(φ₁−φ₄) = 3/10, h = 1/2.
K_κ = s² + κ s[X₁/(h²λ₁) + X₂/(2h⁴yλ₂) + X₃/(6h⁵y²λ₃)] + κ²X₃X₁/(24h⁷y²λ₁λ₃), s = 1 − κ(z₁+z₄), X₁ = τ₁(θ_F1), X₂ = τ₂(u₀) with the prefactor
∏(A: u₀−θ−h; F: u₀−θ), X₃ = reduced τ₃(u₀) (vanishing factor u₀−θ_F1−h removed), u₀ = θ_F1+h, λ_k the product of the Lax constants over the pairs 2 and 3.
Port tests of the fused transfer matrices on vectors [exact/num]: ket and bra, k = 1,2,3 vs explicit matrices; **X̂₃X₁ = 24h⁷λ₁λ₃P exactly at L = 2**; K₁ on the six L = 3 data
[1,0,1] HW states = 3.65685, −7.65685, 14.6580 (×2), 1.842 (×2) (as in the data) [num].

**Procedure** (notebook cells 32–37, 70-digit fixed precision, ε = 1/100, 1/200, 1/400, 1/800; four cells per ε: exact 141×141 sector matrix of t₁^Λ(3/11) by vector propagation,
machine eigensystem + selection by projection on the six data states + inverse iteration, exact c_k (27 words, exact solve at that ε) and bras for τ's, assembly).
Sector (1,0,0,−1), n = 2, 141 states; only the six states selected by projection on the data [1,0,1] HW states are assembled (projection ≥ 0.9998 at ε = 1/100).
C_direct = K·⟨W_a^3|Ψ⟩/(1−κz₁)² (Ψ unit norm, K from the Q-data), C_Q = K·(a/κ₀)²(λ₁u₄)^{−3}⟨0|Ψ⟩Σc_kτ^k/(1−κz₁)².

| data state | |C_direct| at ε = 1/100, 1/200, 1/400, 1/800 | extrapolated (cubic in ε) |C_direct| | extrapolated |C_Q| | |C| of the data |
|---|---|---|---|---|
| #17 | 0.167762, 0.169444, 0.170285, 0.170706 | 0.17112641 | 0.17112641 | 0.17112641 (C = −0.171126) |
| #18 | 0.0261604, 0.0266461, 0.0268912, 0.0270143 | 0.027137747 | 0.027137747 | 0.027137747 (−0.027138) |
| #19, #20 (conj. pair) | 0.0441177, 0.0450046, 0.0454515, 0.0456759 | 0.045900846 | 0.045900846 | 0.045900846 (−0.035872 ∓ 0.028637i) |
| #21, #22 (conj. pair) | 0.151340, 0.152414, 0.152950, 0.153217 | 0.15348408 | 0.15348408 | 0.15348408 (0.129924 ∓ 0.081713i) |

* **C_Q reproduces C_direct to relative 1e-41 … 1e-44 at every ε and every state** [hp] (the F-only 27-word form carries the full ⟨W_a^3|Ψ⟩; the ε⁻⁹ cancellations in ⟨0|Ψ⟩ are
  absorbed at 70 digits).
* The ε → 0 limits of both C_direct and C_Q agree with the six data values to the 8 printed digits (|C|; the phase of Ψ is arbitrary, the conjugate pairs have equal moduli).
  The convergence is linear-plus-higher in ε (|C| − |C_data| ≈ −3.4e-3 at ε = 1/100 for #17, halving with ε), as at L = 2.
* The two L = 3 states with complex conjugate τ's (#19/#20, #21/#22) are not degenerate for t₁^Λ(u) at finite twist (complex conjugate eigenvalues), so the selection is unambiguous.
* [exact] The closed form c_k reproduces the exact solve at two further points: (h, θ, z) = (2/3, (3/11,−5/7,9/13,−13/9,4/9,−2/5), (7/3,3/2,−2/5,5/6)) and the ε-chain point ε = 1/100
  (ranks 27/27, closed form == exact c_k) — notebook cell 37. So the even bridge at L = 3 can be assembled from the closed form without solving a linear system.

Even bridge at L = 3 in Q-data:  C(ε→0) = lim P_κ(τ)·(a/κ₀)²(λ₁u₄)^{−3}·⟨0|Ψ_ε⟩·ũ-form(τ₁(θ_F1),τ₁(θ_F2),τ₁(θ_F3))/(1−κz₁)², with the F-form = (B_{z₂}B_{z₃})(Ŷ)/κ-normalised
words above. As at L = 2 it is Q-complete up to the SoV norm ⟨0|Ψ⟩/‖Ψ‖.

## 5. Failures, caveats, kernel incidents

* The shared kernel was busy for long stretches by the other agent; one of *my* cells (fused-transfer test) hung because `fdBL` lacked the k = 1 key (Missing-key
  message flood) and ran > 6 min; I cancelled only that operation (own operation id) with `wolfbook_cancelOperation(mode = abort)`; the kernel process was not touched. Fixed (k = 1 blocks added, `TimeConstrained` guard).
* `MemoryInUse[]` was 2.4 GB at the start (other agent) and 4.2–4.7 GB later; my own symbols total < 10 MB (`ByteCount`), so the 3 GB limit was exceeded by the shared baseline, not by this work.
* L = 3 only (and L = 2 test). The model ↔ chain equivalence is checked exactly at L = 3 only; the closed form at L ≥ 4 is a model statement.
* The data comparison for the degenerate pairs (19/20, 21/22) is by modulus; selection of the twisted eigenvectors is by projection on the data group (overlap ≥ 0.9998 at ε = 1/100).
