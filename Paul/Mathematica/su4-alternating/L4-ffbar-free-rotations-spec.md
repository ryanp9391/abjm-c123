# Spec: alternating chain f̄ f f̄ f (L = 4) — rotations with free parameters from SoV localisation, à la arXiv:2603.26617 (2026-10-11)

> **Moved 2026-10-11 from gaudin-cg (`Paul/Mathematica/new-conventions/`).** "Session scratchpad" references: scripts are now
> `checks/`, cell dumps `ref/`; paper text extracts were not copied (arXiv 1810.10996, 2002.12341, 2603.26617).
> HANDOFF.md / CLAUDE.md references are to the gaudin-cg repo. Results and status: `README.md` in this folder.

**UPDATE (Paul): companion twist G_c(z), starting covector = the plain SoV vacuum ⟨0| = e₁ᵀ⊗e₄ᵀ⊗e₁ᵀ⊗e₄ᵀ, NO extra global g (Ω = ⟨0|, g = 1); families are polynomials in G_c acting on ⟨0|; generic g only as a cross-check.**

Delegated to a wolfbook-builder (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L4-ffbar-free-rotations.wb` (self-contained, first code cell clears Global`; build/run as before:
JSON, `code` CLI, sources in K1, results as markdown cells; no native outputs). **Finite generic twist only — NO untwisting / ε.**

## 0. Paul's question and the 2603 mechanism

Paul: "look at the 2603 paper; I would like rotations which depend on some free parameters." 2603.26617 (Bargheer–Bercini–
Lefundes–Ryan, App. F, G, M; text extract `rv/bblr2603.txt` in the session scratchpad, App. G from line ~1073, App. M ~1348):
**SoV localisation (G2):** for any global rotation ⟨Z^L|K, the product ∏_{α≤ℓ} t(θ_α) of suitably normalised transfer matrices
acts only on the first ℓ sites, as the transfer matrices of the length-ℓ chain: ⟨Z^L|K ∏ t(θ_α) = ⟨Z^ℓ|K ∏ t^{(ℓ)}(θ_α) ⊗ ⟨Z^{L−ℓ}|K.
Linear combinations then give rotations of the first ℓ sites with a free parameter: (G5) ⟨K^L|(1 − κ t[1,0](θ₁) + κ² t[1,1](θ₁))
= ⟨K'| ⊗ ⟨K^{L−1}| (resolvent torus), the second parameter ω deforms the global vacuum W_ω, and the overlap of the two-rotation
covector ⟨W_{κ,ω} ⊗ W_ω^{L−ℓ}| with any eigenstate is a Q-determinant (G11)–(G18). We want the same for the alternating chain.

**Naming (Paul's, use it in the notebook):** a *fundamental* site is one whose SoV layer is built with t₁ⁿ, n = 0..3; an
*anti-fundamental* site uses t_n, n = 0..3. In the repo's new Lax u + hΣe_ij⊗ρ(E_ij) the anti-fundamental site is
ρ_A(E_ij) = e_ij (HW e₁, group acts by g) and the fundamental site is ρ_F(E_ij) = −e_ji (HW e₄, group acts by g^{−T}) — the
λ ↔ λ* swap of the new conventions. Chain: sites 1..4 = A F A F, T = ℒ⁽¹⁾ℒ⁽²⁾ℒ⁽³⁾ℒ⁽⁴⁾ (site 1 leftmost), twist on the right,
companion G(z), t_a as quantum minors. L = 2 SoV basis (previous notebook `L2-ffbar-SoV-localisation.wb`, cells 18–21):
⟨a,k| = ⟨0|t_a(θ₁−h)t₁(θ₂)^k, ⟨0| = e₁ᵀ⊗e₄ᵀ.

## 1. Main-session findings (numpy, 10⁻¹⁵, to make exact)

Ω = ⟨0|ρ(g) at L = 4 (⟨0| = e₁ᵀ⊗e₄ᵀ⊗e₁ᵀ⊗e₄ᵀ), Ω₂ the same on two sites (so Ω = Ω₂⊗Ω₂). Point: h = 1/2,
θ = (2/7, −3/5, 0.731→ use 8/11, −1.37→ use −11/8), z = (3/2, 5/4, 2/3, 4/7) (exact rationals in the notebook).
(L1) ⟨Ω|t_a(θ₁−h) t₁(θ₂)^k = c_{ak} ⟨Ω₂|t^{(2)}_a(θ₁−h) t^{(2)}₁(θ₂)^k ⊗ ⟨Ω₂| **only for a ∈ {0,3}, k ∈ {0,1}** (4 words); all
other SoV words fail (residual 0.03–0.9). So the (G2) property holds with ONE nontrivial operator per site.
(L2) Singles that localise on block {1,2}: t₃(θ₁−h), t₁(θ₁)⁻¹, t₁(θ₂). Not: t₁(θ₁), t₃(θ₂+3h), t₁(θ₂)⁻¹, t₂ at θ₁, θ₁±h,
θ₂, θ₂+h, θ₂+2h. Localising pairs: t₁(θ₁)⁻¹·t₁(θ₂) (block rotation by G^{−T}, exact, previous notebook cell 33) and the
trivial ones (t₁(θ₁)t₃(θ₁−h), t₁(θ₂)t₃(θ₂+3h) — check whether these are ∝ 1 on Ω, qdet-type).
(L3) Previous notebook: ⟨Ω|t₁(θ₁)⁻¹ = c₁⟨Ω|ρ₁(G^{−T}), c₁ = 1/(hθ₁₂(θ₁₃+h)θ₁₄); at L = 2 the site-2 operator carries a singlet
term: ⟨Ω₂|t₁(θ₂) = hθ₁₂⟨Ω₂|ρ₂-plain(G) − h²⟨Ω₂|G₂|S⟩⟨S|.
**Hence the first free-parameter family:** ⟨Ω|(1 − κ c t₃(θ₁−h)) = e₁ᵀg(1 − κG^{−T}) ⊗ (sites 2..4 rotated by g) — two distinct
rotated vacua, free κ (any linear combination of two site-1 covectors is again e₁ᵀ times a group element, since the HW orbit
of an A site is C⁴∖0), overlap (1 − κ c τ₃(θ₁−h))⟨Ω|Ψ⟩.

## 2. Tasks

**Task 0 — setup.** L = 4 chain as above (256-dim; exact rationals if feasible, else 60-digit numerics with exact rational
inputs), the L = 2 block chain, twisted eigenbasis of the L = 4 chain (simple spectrum check), Ω for an exact rational g.

**Task 1 — localisation inventory (exact).** Make (L1)–(L3) exact; add: fused T_λ (|λ| ≤ 6, ≤ 4 rows) at θ₁ + mh and θ₂ + mh
(m = −5..5) and their inverses — does any of them localise on site 1 or block {1,2} AND act there as a DIFFERENT power of the
twist (G^{−2T}, G^{−3T} on site 1; G², G⁻¹ … on site 2)? That is what would give more than one parameter per site (the full
torus C[G] on a site is 4-dim). Also: the right end (site 4, block {3,4}) by symmetry; one table: operator → localised action
(which power of G on which site, singlet term or not) → constant in closed form in θ, h.

**Task 2 — free-parameter rotation families.** (a) The κ-family above, exact, with the closed-form constant and the overlap
identity on all eigenstates. (b) Block {1,2}: the 4-dim span of the localising words on Ω₂; find all its elements that are
(i) product covectors (rank 1 as 4×4) — i.e. sites 1, 2 rotated independently, (ii) pure block rotations ⟨Ω₂|ρ₁₂(K)| — and
parametrise them (κ₁, κ₂ or a single κ; is the singlet term an obstruction, and can the two parameters cancel it?). Give the
rotations as explicit group elements of the form g·f(κ, G). (c) Both ends at once: ⟨Ω|(1 − κ c t₃(θ₁−h))(1 − κ′ c′ (site-4 operator))
— does it factorise (sites 1 and 4 deformed independently, sites 2, 3 untouched)? (d) The ω-type parameters: global torus
ρ(X), X ∈ C[Gᵀ] invertible, is in the Bethe algebra (commutes with t_a; verify) — express ρ(X) or at least its eigenvalue on Ψ
(expected ∏ x_i^{μ_i} in the diagonal frame) and combine: ⟨Ω|ρ(X)(1 − κ c t₃(θ₁−h))…. Report the total count of free
parameters in the family and which relative rotations between sites are reachable.

**Task 3 — reach of the family (finite twist).** Rank of the span of the covectors in each family (over random parameter values)
inside the 256-dim space and inside each weight space of the untwisted gl(4) action; compare with the global-rotation family
alone (the global family is blind to all non-top irreps; the point of the free parameters is to reach the others — report which
irreps of (1⊕15)⊗(1⊕15) become visible, using the untwisted invariant projectors, NO ε).

**Task 4 (stretch) — Q-form.** For the κ-family the overlap is (1 − κ c τ₃(θ₁−h)) × ⟨Ω|Ψ⟩; express τ₃(θ₁−h) through the Q's
(quantum Wronskian at the point) and, if time permits, ⟨Ω|Ψ⟩ in the L = 4 SoV basis (expansion support, as in 2603 (F3)).

**Summary cell:** table of claims with [exact]/[numeric]/[conjecture] and cell numbers; the explicit free-parameter families
(group elements per site, as functions of κ, ω/X, g, G); parameter count; reach.

## 3. Gotchas
Multi-statement cells need `;`; a comment containing `m*)` closes early; memoised definitions cache bad results;
`wolfbook_evaluateExpression` truncates output (~300 chars) — write long results to files and Read them; `.` binds tighter
than `/`. Don't touch `paper/*.tex`, don't commit.
