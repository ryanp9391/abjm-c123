# Spec: su(4) alternating chain f ⊗ f̄, L = 2 — SoV basis (Ryan–Volin), exact localisation, rotations acting on few sites (2026-10-11)

> **Moved 2026-10-11 from gaudin-cg (`Paul/Mathematica/new-conventions/`).** "Session scratchpad" references: scripts are now
> `checks/`, cell dumps `ref/`; paper text extracts were not copied (arXiv 1810.10996, 2002.12341, 2603.26617).
> HANDOFF.md / CLAUDE.md references are to the gaudin-cg repo. Results and status: `README.md` in this folder.

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-ffbar-SoV-localisation.wb` (self-contained, first code cell clears Global`; build/run
as every recent new notebook: write the .wb JSON, open with the `code` CLI, evaluate each code cell's exact saved source in
K1 with `wolfbook_evaluateExpression`, store printed results as markdown cells under each code cell, save; say in the
summary that there are no native outputs).

## 0. Paul's question

"Consider the su(4) alternating spin chain fundamental ⊗ anti-fundamental ⊗ fund ⊗ anti-fund, start with L = 2. Use the
companion-twist-frame construction (Ryan–Volin) to build the SoV basis. Study SoV localisation in this context. Can we
build global rotations acting on only a few sites?"

**References (read the relevant parts; text extracts are in the session scratchpad, `rv/rv2018.txt`, `rv/rv2020.txt`):**
- RV2018 = arXiv:1810.10996 (companion twist frame; rectangular reps (S^A), SAME rep on every site — does NOT cover f ⊗ f̄,
  since f has A = 1 and f̄ has A = N−1 = 3 and the shortening condition (3.39) depends on A).
- RV2020 = arXiv:2002.12341 (any finite-dim irrep at EACH site; this is the one to use). Key formulas: B-operator (2.25)
  with the (modified) companion twist (2.20), GT patterns / dual diagonals (3.13)–(3.15), SoV basis (4.12)
  ⟨Λ| = ⟨0| ∏_α ∏_{k=1}^{n−1} τ^{n−k−1}[T_{μ̄_k^α}(θ_α + ħ ν^α_{k+1})] (τ^r: T_ij → T_{i+r,j+r}, lower k to the left),
  equivalently (1.2) with ratios of full gl(n) transfer matrices T_{λ̄_k}/T_{(ν_{k+1}^k)}, and the Bäcklund form (1.3);
  wave functions (1.1)/(5.x) as products of Slater determinants of q's. §5 (txt lines ~1174–1700).

**Convention dictionary (verify, don't trust).** RV Lax is u − ħ Σ E_ij ⊗ π(E_ji); ours is u + h Σ e_ij ⊗ ρ(E_ij). With
ħ = h these agree for π = ρ* (π(E_kl) = −ρ(E_lk)), so **our rep λ is RV's rep λ\***: our site-1 f is RV's f̄ and our site-2
f̄ is RV's f. RV's ⟨0| (dual of the lowest weight state of π) is then the dual of the HIGHEST weight vector of ρ — e₁ᵀ for
our f, e₄ᵀ for our f̄ — consistent with the repo's ⟨x₀| = e₁ᵀ. RV order the monodromy ℒ_L…ℒ_1; ours is ℒ⁽¹⁾ℒ⁽²⁾ (site 1
leftmost in the aux product) — irrelevant for the basis (product over α) but keep it in mind for B. The repo's gl(3)
new-convention B (`L2-rotation-sov.wb`, `Bop2`): B(u) = T₁₁(u) qmin_{12,12}(u+h) + T₂₁(u) qmin_{12,13}(u+h), untwisted T —
which is RV (2.27) for N = 3 verbatim. So the gl(4) candidate is
  B(u) = Σ_{k₁∈{1,2,3}} Σ_{K₂⊂{1,2,3},|K₂|=2} T_{k₁,1}(u) · qmin_{K₂; (1, k₁+1)}(u+h) · qmin_{(1,2,3); (1, K₂+1)}(u+2h)
(qmin_{rows; cols}, rows antisymmetrised, column k at v−(k−1)h, as `qminor` in `L2-rotation-sov.wb`; K₂+1 = add 1 to
each entry; ordering of the column list as written). Determine the evaluation points of the SoV basis by the B-eigen test,
using the dictionary only as the starting guess.

## 1. Setup and conventions

- Sites: site 1 = f, ρ_f(E_ij) = e_ij (4×4 matrix unit), HW e₁. Site 2 = f̄, ρ_f̄(E_ij) = −e_ji, HW e₄ (weight
  (0,0,0,−1)). Group: ρ_f(g) = g, ρ_f̄(g) = g^{−T}. Optional cross-check: f̄ ≅ Λ³ = [1,1,1,0] (`JJ` GT code) up to a scalar
  shift of θ₂ — one cell, do not use it downstream.
- Lax ℒ⁽ˢ⁾(u) = Σ_ij e_ij ⊗ ((u−θ_s)δ_ij + h ρ_s(E_ij)); T(u) = ℒ⁽¹⁾ℒ⁽²⁾ (aux product, site 1 left); twist on the right
  T^G = T·G; t_a(u) = Σ_{|I|=a} qmin_{I,I}[T^G](u), a = 1..4; fused T_λ by CBR as in the repo.
- Companion twist G(z) = {{c₁,−c₂,c₃,−c₄},{1,0,0,0},{0,1,0,0},{0,0,1,0}}, c_k elementary symmetric (as `gcMat` in
  `L2-gl4-rect-untwist.wb`); P_aux(z) and the diagonal frame as there.
- Covectors are rows acted on from the right. K = Σ_ij e_ij ⊗ e_ij on f ⊗ f̄ (= |S⟩⟨S| with S = Σ_i e_i ⊗ e_i, the
  invariant; check ⟨S|ρ_f(g)⊗ρ_f̄(g) = ⟨S|), P = Σ e_ij ⊗ e_ji.
- Point P1: h = 1/2, θ = (2/7, −3/5), z = (3/2, 5/4, 2/3, 4/7). Exact rationals throughout (16-dim space; the L = 4
  preview is 256-dim, still fine exactly or at high precision). Second point for any closed form you claim.
- Reference code: `ref-gl4-rect-untwist.txt` (gl(4) companion matrix, P_aux, quantum minors `qminL`, `tOne`) and
  `ref-rotation-sov.txt` (gl(3) L = 2: `qminor`, `tAt`, `Bop2`, SoV basis, B checks) in the session scratchpad.

## 2. Main-session findings to verify first (numpy, 10⁻¹⁵; NOT yet exact)

(F1) **qKZ form at the inhomogeneities, L = 2.** t₁(θ₁) = h·G₁ᵀ·(θ₁₂ − hK₁₂) and t₁(θ₂) = h·G₂·(θ₁₂ − hK₁₂), where G₁ᵀ is
the plain matrix Gᵀ on site 1 (= ρ_f(Gᵀ)) and G₂ the plain matrix G on site 2 (= ρ_f̄(G^{−T})). Reason: ℒ_f(0) = hK_{a1} is
rank 1 and ℒ_f̄(0) = −hP_{a2}, so the aux trace transports the twist onto one site.
(F2) **L = 2 spectrum is almost free.** 12 of the 16 twisted eigenvectors are rank-1 as 4×4 matrices (product states,
expected e_i ⊗ ē_j, i ≠ j, in the diagonal frame); the other 4 are rank 4 (zero-weight states mixed with the singlet).
(F3) **L = 4 alternating chain f f̄ f f̄, qKZ form** (fit over orderings; first match reported, uniqueness not checked):
  t₁(θ₁) = h G₁ᵀ R₁₄ R₁₃ R₁₂, t₁(θ₂) = −h R₂₃ R₂₄ G₂ R₂₁, t₁(θ₃) = h R₃₂ R₃₁ G₃ᵀ R₃₄, t₁(θ₄) = −h G₄ R₄₁ R₄₂ R₄₃,
  R₁₂ = θ₁₂ − hK, R₁₃ = θ₁₃ + hP, R₁₄ = θ₁₄ − hK; R₂₁ = θ₂₁ + hK, R₂₃ = θ₂₃ + hK, R₂₄ = θ₂₄ − hP; R₃₁ = θ₃₁ + hP,
  R₃₂ = θ₃₂ − hK, R₃₄ = θ₃₄ − hK; R₄₁ = θ₄₁ + hK, R₄₂ = θ₄₂ − hP, R₄₃ = θ₄₃ + hK (K, P on the indicated pair).
  f–f̄ R-matrices are θ ± hK: scalar except in that pair's singlet channel; f–f and f̄–f̄ are θ ± hP.
(F4) **End-site localisation on any globally rotated vacuum, L = 4.** For Ω = ⟨0|ρ(g), ⟨0| = e₁ᵀ⊗e₄ᵀ⊗e₁ᵀ⊗e₄ᵀ, generic g:
  ⟨Ω|t₁(θ₁)⁻¹ ∝ ⟨Ω|ρ₁(G^{−T}) and ⟨Ω|t₁(θ₄)⁻¹ ∝ ⟨Ω|ρ₄(G⁻¹) (residuals 10⁻¹⁵), because every R acts on Ω as a scalar (Ω is
  in the adjoint channel of each f–f̄ pair and symmetric in each f–f / f̄–f̄ pair). Interior sites 2, 3: no (residual ~0.8).
  Consequence: ⟨Ω ρ₁(G^{−T})|Ψ⟩ = (scalar)·⟨Ω|Ψ⟩/τ₁(θ₁) on every eigenstate — a single-site rotation, composed with an
  arbitrary global one, is a Q-number. At L = 2 both sites are end sites. Operator-level: the block rotation
  ρ₁(Gᵀ)ρ₂(Gᵀ) on sites 1,2 of the L = 4 chain does not commute with t₁(u) (so it is not in the Bethe algebra) and
  t₁(θ₁)t₁(θ₂)⁻¹ ≠ it; at L = 2, t₁(θ₁)t₁(θ₂)⁻¹ = the global twist rotation exactly (expected from F1).

## 3. Tasks (in order; one result markdown cell under each code cell; summary cell at the end)

**Task 0 — setup.** Reps, invariance of S, Lax, T, t_a (a = 1..4) at P1; quantum determinant t₄ (closed form);
commutativity [t_a(u), t_b(v)] = 0; twisted eigenbasis (16 states, simple spectrum?), weights in the diagonal frame;
decomposition 1 ⊕ 15 of the untwisted chain and the projector onto the singlet.

**Task 1 — exact localisation (F1, F2).** (a) Prove F1 exactly, symbolic in h, θ₁, θ₂, z (or exact at P1 + a second point if
symbolic is slow). (b) Scan: for which a = 1..4, m = −4..4 is t_a(θ_α + mh) of qKZ/local form (twist-on-one-site ×
R-matrices), and the same for the fused T_λ with |λ| ≤ 4 at θ_α + mh — report a table. (c) F2 exactly: closed forms of the 12
product eigenstates and their t₁ eigenvalues, and of the 4 mixed ones (a 4×4 problem in the zero-weight + singlet space).

**Task 2 — gl(4) B-operator.** Implement the candidate B(u) of §0 with untwisted T (also the MCT variant (2.25) with
w's, as a check of z-independence/w-dependence). Checks: polynomial degree; [B(u), B(v)] = 0; diagonalisable; eigenvalues
∏(u − x) with x on θ_α + hℤ; spectrum simple or not (RV2020 fn. 3: non-degenerate for f/f̄-type reps); ⟨0| = e₁ᵀ⊗e₄ᵀ
is a left eigenvector. Compare the spectrum with GT dual diagonals of each site (via the λ ↔ λ* dictionary). If the candidate
fails, say exactly how (and try the transposed/reversed-index variant before anything else).

**Task 3 — SoV basis (RV2020 (4.12)/(1.2)).** Build the 16 covectors ⟨Λ| = ⟨0| × (transfer matrices at θ_α + shifts); fix the
shifts and the τ-embedded sub-chain transfer matrices by the B-eigen test (each ⟨Λ| must be a left B-eigenvector with the
predicted roots). Checks: z-independence, linear independence (rank 16), triangularity relative to the GT basis.
Wave functions: ⟨Λ|Ψ⟩/⟨0|Ψ⟩ = product of transfer-matrix eigenvalues (automatic); then, as a stretch, the Slater-determinant
Q-form (1.1) with gl(4) Baxter Q's (dual Q's of the repo conventions; Baxter equation from the t_a). Report which separated
variables are dynamical per site (f vs f̄ sites differ: RV-f has three single-box layers, RV-f̄ one column layer).

**Task 4 — rotations acting on few sites (the core question).**
(a) F4 at L = 2 exactly: for Ω = ⟨0|ρ(g) (generic exact rational g), ⟨Ω|t₁(θ₁)⁻¹ = c₁⟨Ω|ρ₁(G^{−T})?, ⟨Ω|t₁(θ₂)⁻¹ =
c₂⟨Ω|ρ₂(G⁻¹)?, c's in closed form (expected 1/(h θ₁₂)); verify the overlap identity on all 16 eigenstates.
(b) Reach: which single-site covectors ⟨Ω|ρ₁(X)| (X ∈ GL(4)) are obtainable from ⟨Ω| by words in t_a(θ_α + mh)^{±1} (all a, m
from Task 1(b)), exactly or up to singlet corrections (⟨S|-terms)? In particular: powers/polynomials of the twist on one site
(the torus commuting with G), and what the singlet correction is once the rotated covector is no longer in the adjoint channel.
(c) Generic site-1 rotation ⟨e₁ᵀρ(g₁) ⊗ e₄ᵀ| and generic two-site ⟨e₁ᵀρ(g₁) ⊗ e₄ᵀρ_f̄(g₂)|: number of SoV basis covectors in
their expansion (support), versus the twist-torus ones of (a)/(b); and the rank of the site-1-rotation family over g₁.
(d) CANCELLED (Paul, 2026-10-11): no untwisting / small-ε work in this thread. Everything stays at finite generic twist.

**Task 5 — L = 4 preview (f f̄ f f̄, 256-dim).** Verify F3 exactly (or at 50 digits) and check alternative orderings; F4 for
end sites; for interior sites find what ⟨Ω|t₁(θ₂)⁻¹ is (expected ⟨Ω|ρ₂(G⁻¹)·(R-dressing to sites 3,4)) and whether a word in
t₁(θ_α)^{±1} (length ≤ 3) maps ⟨Ω| to a covector rotated only on a block of sites ({1,2}, {3,4}, {2,3}) — exactly, or modulo
singlet-channel terms. Report the obstruction precisely if none does.

**Summary cell:** a table of every claim with status [exact]/[numeric]/[conjecture] and the cell number; the answer to
"can we build global rotations acting on only a few sites?" in one paragraph; any structural surprise flagged.

## 4. Gotchas

- Multi-statement cells need `;`; a comment containing `m*)` closes early and the cell silently returns Null.
- Memoised definitions cache bad results — re-run downstream after a fix.
- Do not reuse the shared-kernel names of other notebooks without the ClearAll in cell 1.
- Don't touch `paper/*.tex`, don't commit.
