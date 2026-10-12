# ABJM three-point functions from localising transfer matrices — status 2026-10-11

Goal (Paul): combine the ABJM note (`StructuringABJM.pdf`, repo root), Bargheer–Bercini–Lefundes–Ryan
arXiv:2603.26617 ("2603") and our rotation experiments into a formula for the tree-level structure constant of the
note's §3 written **in terms of transfer matrices** — i.e. C = K(τ)·⟨W^L|Ψ⟩ as in 2603 — ignoring Q-functions for now.

**Bottom line so far (L = 2 pairs = 4 sites, numerics at one generic point):**

1. The 2603 mechanism transfers verbatim: reference = an integrable global rotation ⟨W^L| of the vacuum, O₂ = a
   twist rotation 1 − κΛ of the first ℓ fields of the reference, produced by a polynomial K in transfer matrices.
2. **Odd bridge, one field (ℓ₁₂ = 1): works exactly**, C ∝ (1 − z₄/(c₁τ₁(θ₁)))/(1 − z₄/z₁)·⟨W₋₁^L|Ψ⟩ (§5.5),
   verified on on-shell highest-weight states through the untwisting limit (§5.7).
3. **Even bridge (the note's (3.6), one pair): not reached by localising words.** On a pair block the localising words
   give only three rotation types; the missing one is Λ acting on the F site alone (§5.3).
4. Paul's highest-weight trick (global Borel rotations are free for HW Ψ) makes (3.6) **exactly** reachable at finite
   twist with a single word, but the required rotation parameter b = z₁/(z₁−z₄) diverges as z → 1 and the untwisting
   limit gives the wrong, state-dependent answer (§5.6). 2603 avoids this with a continuous κ-family (κ → 1 together
   with z → 1); for ABJM pairs that family is exactly what is missing.

Evidence tags as in `README.md`: **[numpy, 10⁻ⁿ]** = double-precision check in `checks/`; nothing here is in a
notebook yet. Point: h = 1/2, θ = (2/7, −3/5, 8/11, −11/8), z = (3/2, 5/4, 2/3, 4/7) unless stated.

---

## 1. Sources and what was taken from them

### 1.1 The ABJM note (`StructuringABJM.pdf`, draft, 7 pp.)

- (3.3) BPS operators O₂ = tr(Y₁†Y⁴)^L, O₃ = tr(Ỹ₁†Ỹ⁴)^L with Ỹ₁† = Y₁† − Y₄†, Ỹ⁴ = Y¹ + Y⁴.
- (3.4)–(3.5) O₁ is an excitation of the vacuum tr(Y¹Y₄†)^L (2L fields).
- (3.6) C^{•∘∘}/C^{∘∘∘} = 𝒩₁₂₃ (⟨Y₁†Y⁴|^{L−M} ⊗ ⟨Ỹ₁†Ỹ⁴|^{⊗M}) |Ψ⟩: the first L−M pairs of O₁ contract with O₂, the
  last M pairs with O₃.
- (4.1) Y^a ↦ e^a, Y_a† ↦ e_a. (4.8) integrable boundary state ⟨W|(τ̂₄ − τ̂₄̄) = 0. (4.9) ⟨W_a| = (e₁ + a e₄) ⊗ (e⁴ − a e¹).

### 1.2 arXiv:2603.26617 (LaTeX source read from arXiv; key places)

- Chain of 6's (∧²C⁴), Lax ℒ_jk(u) = uδ_jk − iE_kj, **diagonal twist** T_jk(u) = z_j Σ ℒ⁽¹⁾…ℒ⁽ᴸ⁾ (App. B).
- §II: O₂ = Tr 𝒲_{ω,κ}^{L₂}, O₃ = Tr 𝒲_ω^{L₃}, each an SU(4) rotation of ⟨Z| (eq. `WAdef`);
  C_ℓ = ⟨𝒲_{ω,κ}^ℓ ⊗ 𝒲_ω^{L−ℓ}|Ψ⟩/‖Ψ‖ (eq. `twistedOverlap2`).
- App. A: polarisations p^{ω,κ}_{ij} = p^ω_{ij}(1−κz_j)/(1−κz_i); ν_j are ratios of twist differences, e.g.
  ν₁ = (z₁₃z₂₄/(z₁₂z₃₄))^{1/2}.
- §V "Structure Constant": localisation on any global rotation of the vacuum,
  ⟨𝒲_ω^L|∏_{k≤ℓ} t̂_{s_k}(θ_k) = ⟨𝒲_ω^ℓ|∏ t̂^{(ℓ)}_{s_k}(θ_k) ⊗ ⟨𝒲_ω^{L−ℓ}| (eq. `localized`), s_k ∈ {+,0,−};
  ⟨𝒲_{ω,κ}^ℓ| = ⟨𝒲_ω^ℓ|K_κ with K_κ a polynomial in t₊, t₀, t₋, so the overlap is K_κ(τ)⟨𝒲_ω^L|Ψ⟩, and K_κ(τ) is
  an su(3) SoV determinant (eq. `Kkappa`).
- §VII: untwisting limit z_i → 1, κ → 1, ω → 1 taken **together**.

**Reading of 2603 used below.** Since conjugating ⟨Z|exp(Σ p E_ij) by D = diag(d_i) multiplies p_ij by d_j/d_i, the
relation p^{ω,κ} = p^ω(1−κz_j)/(1−κz_i) says ⟨𝒲_{ω,κ}| ∝ ⟨𝒲_ω|ρ(D_κ), D_κ = 1 − κ·(twist). So O₂ is the reference O₃
state with the **first ℓ sites rotated by 1 − κΛ** — the "common twist rotation of a block" we built in the companion
frame (README, {1,3} common rotation), written in the diagonal frame.

---

## 2. Conventions (this file)

As in `README.md` (Lax with ρ_A(E_ij) = e_ij, ρ_F(E_ij) = −e_ji, site 1 leftmost, twist on the right; bra site vectors
transform as v → v g on A sites and w → w g^{−T} on F sites). Chain A F A F (L = 2 pairs).

- **Diagonal (physical) frame:** twist Λ = diag(z₁, z₂, z₃, z₄). Bra vacuum per pair ⟨vac| = e₁ᵀ ⊗ e₄ᵀ.
- γ := (twist)^{−T}; in the diagonal frame γ = Λ⁻¹, acting as Λ⁻¹ on A sites and as Λ on F sites.
- **ABJM pairs (assumed identification of the note's (e_a, e^a) on one pair with our (A, F) bra vectors):**

  W_a := (e₁ + a e₄) ⊗ (e₄ − a e₁) = ⟨vac|ρ(1 + aE₁₄) (on one pair), W₀ = vacuum pair (O₂), W₋₁ = O₃ pair.

- **Twist rotation of a pair:** ρ_pair(D) with D = diag(d_i) maps W_a → W_{a d₄/d₁} (up to normalisation). In
  particular D_κ = 1 − κΛ gives W_a → W_{a(1−κz₄)/(1−κz₁)}, and κ = 1/z₄ turns the O₃ pair into the vacuum pair.
  So **(3.6) = ⟨W₋₁^L| with its first L−M pairs rotated by D_{1/z₄}**, the exact ABJM analogue of 2603's O₂.

---

## 3. Companion frame ↔ diagonal frame [numpy, 10⁻¹⁵; `checks/chk_frame.py`]

With V the Vandermonde matrix with columns (z_j³, z_j², z_j, 1)ᵀ (G_c V = V Λ),

  t^{G_c}(u) · R = R · t^{Λ}(u),  R = ⊗_sites (V^{−T} on A sites, V on F sites).

Bras map as ⟨x_diag| = ⟨x_comp|R. Consequences:

- Companion site rotations 1 − κG_cᵀ (A) / 1 − κG_c (F) become 1 − κΛ in the diagonal frame; everything we found in the
  companion frame about localisation is frame-covariant.
- The companion SoV vacuum ⟨0|_c is, in the diagonal frame, a global rotation of the vacuum with **all** eigen-components
  non-zero (e₁ is cyclic for G_c — this is what makes it an SoV vacuum). The physical vacuum has three zero components
  per site, so from ⟨0|_c it would need a degree-3 (Lagrange-projector) twist polynomial on every site; localising
  words give degree ≤ 1 per site (≤ 2 for the {1,3} common rotation). **Hence the reference for the 3-pt function is
  taken directly as a global rotation in the diagonal frame (as in 2603), not built from ⟨0|_c.**
- In the diagonal frame the vacuum itself is a transfer-matrix eigenstate, so it is useless as a reference; the
  reference must be a non-trivial global rotation (Paul).

---

## 4. Companion-frame rotation results of this session (input)

(Also summarised in `README.md`, "Globally rotated reference states" and "Rotations of the 2-parameter state".)

1. **{1,3} common rotation is global.** ⟨0|ρ₁(g)ρ₃(g), g = 1−κG_cᵀ, only sees row 1 of g; g → g′ = 1 − κ(c₁E₁₁ + E₁₂) ∝
   exp(−κ′E₁₂), κ′ = κ/(1−κc₁), fixes the F vacuum, so ⟨0|ρ₁(g)ρ₃(g) = ⟨0|ρ(g′) on all four sites. [algebra]
2. **Good-reference criterion.** A global rotation v⊗w⊗v⊗w (v·w = 0) behaves like ⟨0| (all four sites localise, both
   common rotations exist, same coefficients) iff w G vᵀ = 0. [numpy 10⁻¹⁵, `chk_rotref.py`]
3. **Reach.** v = e₁ᵀp(Gᵀ), w = e₄ᵀq(G): global iff deg p + deg q ≤ 2, good reference iff ≤ 1. ⟨0|P₁₃(κ₁)P₁₃(κ₂) and
   ⟨0|P₁₃P₂₄ are global rotations. [numpy, `chk_rotref.py`]
4. **The 3⊗3 word span.** With X = t₁(θ₁), Y = t₁(θ₃): the nine ⟨0|X^aY^b (a, b ≤ 2) are independent and span exactly
   {e_i⊗e₄⊗e_j⊗e₄ : i, j ≤ 3}; F sites never move. Global rotations inside = v⊗e₄⊗v⊗e₄ with v ∈ span(e₁, e₂, e₃), i.e.
   v = e₁ᵀ(1−κ₁Gᵀ)(1−κ₂Gᵀ) = ⟨0|P₁₃(κ₁)P₁₃(κ₂) (nothing new); also (Gᵀ)² on site 1 alone is inside; v₄ ≠ 0 is not.
   [numpy 10⁻¹⁴, `chk_span13_words.py`]
5. **Generators.** These global rotations are exp(βE₁₂ + β′E₁₃)·(E₁₁ scale) in the gl(3) on indices {1,2,3} (the other
   six gl(3) generators annihilate both vacua). With s = κ₁+κ₂, p = κ₁κ₂:
   e₁ᵀ(1−κ₁Gᵀ)(1−κ₂Gᵀ) = (1 − sc₁ + p(c₁²−c₂)) e₁ + (pc₁ − s) e₂ + p e₃. [algebra]
6. **The 2-parameter state Ω₂ = v⊗e₄⊗v⊗e₄ is only a generic global rotation** (w G vᵀ = κ₁κ₂ ≠ 0). Exhaustive scan
   (t₁, t₂, t₃ at θ_s + mh, m = −2…2, and inverses; single words and pairs; triples of t₁(θ_s)^{±1}, t₃(θ_s−h)^{±1}):
   only end-growing words localise, all as fixed twist rotations, identical to a random global rotation.
   [numpy, `scan_*.py`]

   | word on ⟨Ω₂\| | result |
   |---|---|
   | t₁(θ₁)⁻¹ (∝ t₃(θ₁−h)) | ρ₁(γ) |
   | t₁(θ₁)⁻¹ t₁(θ₂) | ρ₁₂(γ) |
   | t₁(θ₁)⁻¹ t₁(θ₂) t₁(θ₃)⁻¹ (t₃(θ₃−h) can replace t₁(θ₃)⁻¹) | ρ₁₂₃(γ) |
   | t₁(θ₁)⁻¹ t₁(θ₂) t₁(θ₃)⁻¹ t₁(θ₄) | ρ(γ), global torus |
   | t₁(θ₄)⁻¹ | ρ₄(γ⁻¹) |
   | t₁(θ₃) t₁(θ₄)⁻¹ | ρ₃₄(γ⁻¹) |
   | t₁(θ₂)⁻¹ t₁(θ₃) t₁(θ₄)⁻¹ | ρ₂₃₄(γ⁻¹) |

   Free parameters only on a single end site; t₁(θ₁)⁻¹t₁(θ₄)⁻¹ entangles sites 1 and 4. Interior single sites and the
   {1,3}/{2,4} common rotations are lost.

**Lesson for the diagonal frame:** on a generic global rotation, localisation grows from the ends of the chain, as in
2603's eq. `localized`; this is all that is available on the ABJM reference ⟨W_a^L| (which is not a "good reference":
w Λ vᵀ = a(z₄ − z₁) ≠ 0).

---

## 5. Diagonal frame, ABJM reference ⟨W₋₁|^{⊗2} (L = 2 pairs)

Scripts: `checks/abjm_pair_*.py`, `abjm_odd_bridge.py`, `abjm_untwist*.py`.

### 5.1 The L = 2 SoV words of sites 1, 2 [numpy; `abjm_pair_base.py`]

Of the 16 words t_a(θ₁−h) t₁(θ₂)^k (a, k = 0..3) acting on ⟨W₋₁|^{⊗2}, exactly {1, t₃(θ₁−h)} × {1, t₁(θ₂)} localise on
pair 1 (sites 3, 4 stay W₋₁); their span ∩ (pair-1 states ⊗ W₋₁) is 4-dimensional and does **not** contain the
rotated reference ρ₁₂(1−κΛ)⟨W| (residual 2–6·10⁻²).

### 5.2 Cyclic span — why "some K exists" is empty [numpy; `abjm_pair_span.py`]

Library: t_k(θ_s + mh), k = 1, 2, 3, s = 1, 2, m = −2…2, plus inverses where invertible. Six of them annihilate ⟨W|:
t₂(θ₁+h), t₃(θ₁+h), t₃(θ₁+2h), t₂(θ₂), t₃(θ₂), t₃(θ₂+h). ⟨W|·{1, singles} has rank 33; adding products of two gives
rank 54, which is the whole cyclic space ⟨W|·(transfer-matrix algebra) (= number of eigenstates with ⟨W²|Ψ⟩ ≠ 0).
The full W_{a′}-family on pair 1 lies in it (10⁻¹⁵), but backward elimination needs all 54 words: the statement
"a polynomial K exists" is just the selection rule, not a 2603-type formula.

### 5.3 Localising words on pair 1 [numpy; `abjm_pair_classify*.py`, `abjm_pair_localspan.py`]

All words with ≤ 2 library factors per site (435 × 406 combinations): product states on pair 1 ⊗ W₋₁ occur only as
ρ₁(Λ^i on A) ρ₂(Λ^j on F) with

| (i, j) | words | example |
|---|---|---|
| (0, 0) | 208 | 1 |
| (−1, 0) | 52 | t₃(θ₁−h) ∝ t₁(θ₁)⁻¹ |
| (−1, 1) | 12 | t₁(θ₁)⁻¹ t₁(θ₂) = ρ₁₂(γ) |

and **all 320 pair-1-localising words, entangled ones included, span only 4 dimensions**: ⟨W|, ρ₁(Λ⁻¹)⟨W|, ρ₁₂(γ)⟨W|,
⟨W|t₁(θ₂). Not in that span: the vacuum pair ⊗ W₋₁ (residual 0.6), the W_{a′} family, the singlet ⟨S|₁₂ ⊗ W₋₁ (0.4),
and the site-2-only rotation ρ₂(Λ)⟨W| (0.2). No (−2, ·) appears: t₁(θ₁)⁻² etc. do not localise.

**Why this is fatal for even bridges.** Twist rotations rescale the four components of a W pair:
ρ₁(D₁)ρ₂(D₂)W_a = (d₁e₁ + a d₄e₄) ⊗ (f₄e₄ − a f₁e₁), with coefficient vector (d₁f₄, d₁f₁, d₄f₄, d₄f₁). The three available
types give (1,1,1,1), (z₁⁻¹, z₁⁻¹, z₄⁻¹, z₄⁻¹), (z₄/z₁, 1, 1, z₁/z₄); the vacuum pair needs (1,0,0,0), and the second
and third components force all coefficients to vanish. A fourth type, Λ on the **F site alone** (vector
(z₄, z₁, z₄, z₁)), would close the system — and is also what a κ-family on the pair needs: ρ_pair(1−κΛ⁻¹) on a W pair
has F-part ∝ 1 − κ(z₁+z₄−Λ)/(z₁z₄), i.e. contains Λ on the F site alone. ⟨W|t₁(θ₂) does give ρ₂(Λ)⟨W| but with a
singlet admixture on (1,2) (its K-term (wΛ)·v = a(z₄−z₁) ≠ 0), and nothing in the scan cancels it.

The mirror (right end) gives (0,0), (0,−1) [t₁(θ₄)⁻¹], (1,−1) [t₁(θ₃)t₁(θ₄)⁻¹] — same deficiency.

### 5.4 Global torus as transfer matrices [numpy, operator identity; `abjm_odd_bridge.py`]

  t₁(θ₁)⁻¹ t₁(θ₂) t₁(θ₃)⁻¹ t₁(θ₄) = c · ρ(γ)  (Λ⁻¹ on A, Λ on F, every site), c = −6.37197 at the point.

So the global torus is an explicit product of transfer matrices. But combining the two ends does not help:
⟨W|P_L P_R and ⟨W|P_L² (P_L = t₁(θ₁)⁻¹t₁(θ₂), P_R = t₁(θ₃)t₁(θ₄)⁻¹) do **not** localise as products.

### 5.5 Odd bridge, one field — the working case [numpy, exact at any twist; `abjm_odd_bridge.py`]

  ⟨W|t₁(θ₁)⁻¹ = c₁ ⟨W|ρ₁(Λ⁻¹),  c₁ = 1/(h θ₁₂ (θ₁₃+h) θ₁₄)

(end-site localisation, holds on any global rotation), hence the 2603 ℓ = 1 κ-family

  ⟨W₋₁^{⊗L}| (1 − κ t₁(θ₁)⁻¹/c₁) = (e₁ + a(1−κ/z₄)/(1−κ/z₁) e₄)_{site 1} ⊗ …,

and at κ = z₄ site 1 becomes (1 − z₄/z₁) e₁ — the O₂ field — while site 2 keeps the O₃ field:

  ⟨e₁ ⊗ (e₄ + e₁) ⊗ W₋₁^{⊗(L−1)}|Ψ⟩ = (1 − z₄/(c₁ τ₁(θ₁))) / (1 − z₄/z₁) · ⟨W₋₁^{⊗L}|Ψ⟩,

  1/c₁ = h ∏_{F sites s} θ₁ₛ ∏_{A sites s ≠ 1} (θ₁ₛ + h)   (general L; L = 2 pairs: hθ₁₂(θ₁₃+h)θ₁₄).

(**Correction 2026-10-11:** an earlier version of this file, of README and of NEXT omitted the factor 1/(1 − z₄/z₁).)

This is the 3-pt function with bridge ℓ₁₂ = 1 field (allowed when L₁ + L₂ − L₃ ≡ 2 mod 4 in fields). By the mirror,
site 4 alone (an F field) has a κ-family via t₁(θ₄)⁻¹ (seen in the Ω₂ scan; not re-checked on W). Exact at finite twist,
so the untwisting limit is continuous. Not done: the homogeneous limit θ → 0 (c₁ is singular — needs the 2603
"properly normalised t̂"), and ℓ = 3 fields (expected to hit the F-site gap of §5.3).

### 5.6 Paul's highest-weight trick and the untwisting limit [numpy; `abjm_untwist*.py`]

**Trick.** For untwisted highest-weight Ψ (raising generators J_ij, i < j, annihilate the ket; the ket vacuum is HW in
these conventions), any global upper-unitriangular ρ(n) acts trivially: ⟨χ|ρ(n)|Ψ⟩ = ⟨χ|Ψ⟩. With n = e^{bE₁₄},
W_a → W_{a+b}, so

  ⟨vac ⊗ W₋₁|Ψ⟩ = ⟨W_b ⊗ W_{b−1}|Ψ⟩  for all b.

Verified on the full 10-dim HW subspace of the untwisted inhomogeneous chain (b = 0.7, −2.3, 5: 10⁻¹⁶–10⁻¹⁴).
(General version: the pair-1 relation must be a twist rotation D with nDn⁻¹ = (1+E₁₄)s, s in the vacuum-pair
stabiliser; for D = γ this forces n = e^{bE₁₄} with the b below — the Borel freedom gives nothing beyond b.)

**Exact finite-twist formula.** ρ₁₂(γ)W_{b−1} ∝ W_{(b−1)z₁/z₄}, so with **b = z₁/(z₁ − z₄)**

  ⟨W_{b−1}^{⊗2}| t₁(θ₁)⁻¹ t₁(θ₂) = c ⟨W_b ⊗ W_{b−1}|   (residual 10⁻¹⁶ at every ε below),

i.e. F := c⁻¹ τ₁(θ₂)/τ₁(θ₁) · ⟨W_{b−1}^{⊗2}|Ψ⟩ = ⟨W_b ⊗ W_{b−1}|Ψ⟩ exactly, with a state-independent c (extracted
numerically; presumably the README end-block constant c₁₂ times z₄/z₁ — not checked).

**Untwisting fails.** z = 1 + εφ, φ = (0.9, 0.4, −0.3, −1.1) and (1.3, −0.2, 0.5, −0.7), ε = 10⁻³, 10⁻⁴; twisted
eigenstates with ‖J_{i<j}Ψ‖ < 30ε and T := ⟨vac ⊗ W₋₁|Ψ⟩ ≠ 0 (7 states). F/T converges, independently of φ, to
**state-dependent** values: 1, 0.846, 0.360, 0.0462, 0.027 ± 0.335i, and one state 0.23–0.27 (not converged).
Mechanism: b ~ 1/ε multiplies the O(ε) non-HW part of the twisted eigenvector, F − T = Σ_{k≥1} b^k f_k(ε) with finite
b^k f_k.

**Comparison with 2603.** There the O₂/O₃ polarisation ratio is (1−κz_j)/(1−κz_i) and κ → 1 with z → 1, so the analogue
of b is (1−κz₄)/(κ(z₁−z₄)) → (1−z₄)/(z₁−z₄): a ratio of twist differences, finite (direction-dependent, which is
harmless by b-independence), and the limit is continuous. Discrete localising words only give ratios (z₁/z₄)^n → 1,
forcing b → ∞. **So the continuous κ-family on a pair is genuinely needed — the F-site gap of §5.3 is the obstruction.**

### 5.7 Explicit odd-bridge computations, on-shell HW states [numpy; `abjm_bridges.py`, `abjm_bridges_t3.py`]

> **2026-10-12 — the θ = 0 set-up below is the wrong chain (`abjm-bridges.wb` §3).** The ABJM (Minahan–Zarembo) point in
> these conventions is **θ_A − θ_F = 2h**: only there does H_MZ commute with t₁, t̄₃ (10⁻¹⁵), and t₁(θ_F) t̄₃(θ_A−h) = h⁴·S
> (pair shift). At θ_A = θ_F this product vanishes, S leaves the t-algebra, and the t-eigenstates are not H_MZ eigenstates.
> Corrected L = 2 table (t₁, t̄₃ alone resolve all states; all H_MZ eigenstates):
>
> | state | E_MZ | nonzero C_ℓ |
> |---|---|---|
> | [2,0,2] vacuum | 0 | all 1 |
> | [0,2,0] | 16 | — |
> | [1,0,1] parity pair | 12, 12 | C₁ = C₃ = (1 ∓ i√3)/√48, \|C₁\| = 1/√12 |
> | singlet | 20 | C₂ = −1/√40 |
> | singlet | 4 | C₂ = 1/√24 |
>
> Same odd/even selection rule. Odd-bridge formula through the untwisting limit at the ABJM point (θ = θ_ABJM + δθ⁰,
> no shift operator): exact to 10⁻¹¹–10⁻¹⁴, |LHS| → 1/√12. The θ = 0 table, caveat and "Resolved" bullet below are
> kept for the record only.

Set-up (`abjm_bridges.py`, transfer matrices applied to vectors, checked against the matrix code): homogeneous
(θ = 0), untwisted chain; on-shell states = joint eigenvectors of t₁ (and t₃) inside highest-weight ∩ cyclic
(translation by one pair = 1, the trace) of each dominant weight sector. Bra for bridge ℓ (fields): first ℓ sites carry
O₂ fields (A: e₁, F: e₄), the rest O₃ fields (A: e₁−e₄, F: e₄+e₁); ℓ = 0 is ⟨W₋₁^L|, ℓ = 2L the vacuum; with this
normalisation C∘∘∘ = 1. C_ℓ = ⟨B_ℓ|Ψ⟩/‖Ψ‖ (Hermitian norm; overall phase of Ψ free).

**L = 2 pairs: 6 cyclic HW states.**

| Dynkin [p,q,r] | C₀ | C₁ (odd) | C₂ (= note's (3.6), M = 1) | C₃ (odd) | C₄ |
|---|---|---|---|---|---|
| [2,0,2] vacuum | 1 | 1 | 1 | 1 | 1 |
| [0,2,0] | 0 | 0 | 0 | 0 | 0 |
| [1,0,1], t₃-eigenstate 1 | 0 | 0.353553 = 1/√8 | 0 | 1/√8 | 0 |
| [1,0,1], t₃-eigenstate 2 | 0 | 0.204124 = 1/√24 | 0 | 1/√24 | 0 |
| [0,0,0] (a) | 0 | 0 | −0.25 | 0 | 0 |
| [0,0,0] (b) | 0 | 0 | 0.258199 = 1/√15 | 0 | 0 |

- **Selection rule (observed):** the adjoint states are seen **only by odd bridges**, the singlets only by even bridges.
  So odd bridges are not a curiosity: at L = 2 they are the only way to see the [1,0,1] states.
- The two [1,0,1] states are degenerate in t₁ (same τ₁(u)) and are separated by t₃ (eigenvalues 3.0667, 2.3419 for
  t₃(0.41) + 0.3 t₃(−0.77)). `abjm_bridges.py` diagonalises with t₁ only, so its printed [1,0,1] C-values (0.3401,
  0.3016) are in an arbitrary basis — use `abjm_bridges_t3.py`. Closed forms above are numerical recognitions.

**Formula through the untwisting limit.** For each HW state with C₁ ≠ 0, deform to z = 1 + εφ, θ = δθ⁰
(φ = (0.9, 0.4, −0.3, −1.1), θ⁰ = (2/7, −3/5, 8/11, −11/8)), ε = δ = 10⁻², 10⁻³, 10⁻⁴, pick the twisted
eigenstate in the same weight sector with largest overlap, and compare LHS ⟨B₁|Ψ⟩ with the RHS of §5.5:

- LHS = RHS to 10⁻¹²–10⁻¹⁶ at every (ε, δ) — exact, including the general-L c₁ above.
- Untwisting: ⟨W₋₁^L|Ψ⟩ → 0 like ε (6.1·10⁻³, 6.2·10⁻⁴, 6.2·10⁻⁵ for the adjoint) while the bracket
  (1 − z₄/(c₁τ₁))/(1 − z₄/z₁) grows like 1/ε (51, 501, 5000): the product is finite. For the vacuum both are 1.
- **Caveat:** for the degenerate adjoint pair the overlap-matching picked the same twisted state for both, and its LHS
  tends to 0.3121, which is neither 1/√8 nor 1/√24. The matching must use t₃ as well (the twisted limit should
  select joint t₁, t₃ eigenstates); unresolved.
- **[Superseded — artefact of the singular point θ_A = θ_F, see box above]** t₃ alone is not enough — the 1/√8 state is degenerate in t₁, t₂, t₃ with an
  odd-momentum HW [1,0,1] state; only the pair translation separates them, and θ = δθ⁰ breaks it. With pair-periodic θ
  and the twisted shift Ŝ = S∘ρ_{3..L}(Λ) (commutes with t_a at finite twist; plain S does not), the matched states
  converge (overlap → 1), LHS → 1/√8 and 1/√24, formula exact to 10⁻¹¹–10⁻¹³.

**Not done (stopped here):** L = 3 pairs. Needs t₃ in vector form for the degeneracies; the naive anti-fundamental
auxiliary Lax u ∓ h Σ e_ij ⊗ ρ(E_ji) does **not** commute with t₁ or t₃ (`abjm_dual_lax_try.py`), so either fused t₃
applied to vectors (quantum minors, ~270 monodromy passes per vector) or the correct dual Lax (shift/sign to be found).

**Update 2026-10-12 — dual Lax found (`dual-lax-t3.wb`, exact):** t₃(u) = ∏_s g_s(u) · tr[ℒ̄₁(u)⋯ℒ̄_L(u) · det G · G^{−T}],
ℒ̄_A = (u−θ+h) − hP, g_A = (u−θ−h)(u−θ−2h); ℒ̄_F = (u−θ−3h) + hK, g_F = (u−θ)(u−θ−h) (sign gauge (−1)^k on the dual index).
The naive try missed the relative shift 4h between A and F sites. Vector form: one pass, 0.16 s per application at
L = 3 pairs, [t₁, t₃] = 0 there to 10⁻¹⁶. L = 3 pairs is now unblocked. Same for t₂ (`dual-lax-t3.wb` §4, exact at
L = 2, 4): ∧² Lax (u−θ̃_s) + h Σ D(e_ab)⊗ρ_s(E_ab) on the 6-dim aux space (D = derivation action), θ̃_A = θ, θ̃_F = θ + h,
g_A = u−θ−h, g_F = u−θ, twist ∧²G; vector form `t2V` in `abjm-bridges.wb`.

---

## 6. Open problems / next steps

1. **Find the pair κ-family**, i.e. an F-site-only twist rotation on the left end (or a combination cancelling the
   singlet admixture of t₁(θ₂)): other shifts / fused T_λ in other reps / longer words; and **L = 3 pairs** (4096-dim),
   in case two pairs is too short for a block plus reservoir.
2. If found: K_κ for the first pair, then for ℓ pairs; check the untwisting limit with κ → 1 as in 2603.
3. Odd bridges (§5.5, §5.7): fix the degenerate-state matching (use t₃), then L = 3 pairs (all odd ℓ, check the
   selection rule and the formula); ℓ = 3 fields; recognise closed forms; homogeneous-limit form of the formula
   (c₁ singular at θ = 0, only the ratio with τ₁ is finite).
4. Is ⟨W_a^L| integrable in the sense of the note's (4.8) for our twisted chain (which twist constraints, analogue of
   2603's z₁z₄ = z₂z₃ = 1, ω = 1)? Not tested.
5. The overlap ⟨W^L|Ψ⟩ itself (2603's 𝕎_ω determinant) — deliberately left for the Q-function stage.

---

## 7. Script index (`checks/`, run `python checks/<name>.py` from anywhere; scripts `exec` siblings)

| Script | Content | Runtime |
|---|---|---|
| `rotbase.py` | shared setup (chain, t's, `prod`, `local`, `report`); its own run part reproduces `chk_rotref.py`'s first half | ~10 s |
| `chk_rotref.py` | §4 items 1–3 (criterion, iterated common rotations) | ~1 min |
| `chk_span13_words.py` | §4 item 4 (3⊗3 span of X^aY^b) | ~10 s |
| `scan_ops_build.py` → `scan_blocks.py` → `scan_classify.py` → `scan_triples.py`, `scan_ends_identify.py` | §4 item 6 (Ω₂ localisation scan; `scan_ops_build.py` must run first, caches operators in the temp dir) | minutes |
| `chk_frame.py` | §3 frame map | ~5 s |
| `abjm_pair_base.py` | §5.1 | ~10 s |
| `abjm_pair_span.py` | §5.2 | ~1 min |
| `abjm_pair_classify.py`, `abjm_pair_classify2.py`, `abjm_pair_localspan.py` | §5.3 | ~1 min each |
| `abjm_odd_bridge.py` | §5.4, §5.5 | ~10 s |
| `abjm_untwist.py`, `abjm_untwist_ratios.py` | §5.6 | ~1 min each |
| `abjm_bridges.py [Lpairs]` | §5.7: HW cyclic states, C_ℓ for all bridges, odd-bridge formula through untwisting (t₁ only) | ~5 s (L = 2) |
| `abjm_bridges_t3.py` | §5.7: degenerate [1,0,1] pair resolved by t₃ | ~1 min |
| `abjm_dual_lax_try.py` | §5.7: naive 4̄-auxiliary Lax does not commute (negative) | ~1 min |

All scripts were smoke-tested from an unrelated working directory (exit 0) on 2026-10-11.
