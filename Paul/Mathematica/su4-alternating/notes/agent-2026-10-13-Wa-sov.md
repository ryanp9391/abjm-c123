# Agent 2026-10-13: SoV expansion of ⟨W_a^L| (L = 2 pairs) and the twisted-vacuum overlap as a finite Q-sum

Notebook: `Wa-sov-expansion-L4.wb` (symbols prefixed `ws`). Python cross-check / fallback: `checks/session-2026-10-13/wasov/`
(`wa_sov_core.py`, `run1`–`run11.py`, logs; exact Fractions + mpmath 50 digits; machinery of `checks/session-2026-10-12/sovq/sov_l4.py`).
Conventions exactly as in README / `notes/w-a-reference-status.md`: chain A F A F, bras row vectors, T = L₁L₂L₃L₄, twist on the right,
bra words ⟨a₁,k₂,a₃,k₄| = ⟨0|t_{a₁}(θ₁−h)t₁(θ₂)^{k₂}t_{a₃}(θ₃−h)t₁(θ₄)^{k₄}, ket words on |0̄⟩ = e₄e₁e₄e₁, M = Gram⁻¹,
standard point h = 1/2, θ = (2/7, −3/5, 8/11, −11/8), z = (3/2, 5/4, 2/3, 4/7), w = (7/5, 9/8, 3/7, 5/9), w₂ = (11/7, 5/3, 2/9, 13/5).
Tags: [exact] exact rationals/symbolic, [hp] 50–60 digits, [num] double.

**Run status.** Notebook cells 1–14 were first run in the shared kernel; that kernel was then wiped/restarted by another agent (~06:50 UTC), so cells 15–28 were first obtained with the python cross-check in `checks/session-2026-10-13/wasov/` (exact Fractions + mpmath 50 digits, same conventions; python was used as fallback while the kernel was down). After the coordinator restarted the kernel (08:03 local) the whole notebook was re-run from the top and **every number below was reproduced in the notebook** (ranks, supports, dual constants, W_w, sectors, eigenbra property, Task 5 counts, candidate tests); where Mathematica (60 digits) and python (50 digits) differ in precision the Mathematica figure is quoted. The notebook is saved with outputs (except where a cell prints only a Grid preview).

> **Update (follow-up section at the end):** on the weight-(μ,0,0,−μ) states ⟨W_a^L|Ψ⟩ has a closed Q-form at L = 2 (bidegree (2,2) polynomial in the two F-site single-Q ratios, 9 F-only SoV words), which supersedes the "no closed Q-form" statements of §3 and §6; the dense 256-word expansion below is the unique representation on all 256 states.

## 1. Companion-frame image of ⟨W_a^L| [exact, symbolic a]

Per pair (diagonal frame) (e₁ + a e₄) ⊗ (e₄ − a e₁), ⟨B_Λ| ↦ ⟨B_Λ|R⁻¹, R⁻¹ = V^T (A) ⊗ V⁻¹ (F):

  ⟨W_a^L| R⁻¹ = [ (v(z₁) + a v(z₄)) ⊗ (ũ₄ − a ũ₁) ]^{⊗L},  v(x) = (x³,x²,x,1), ũ_j = row j of V⁻¹,

coefficient exactly 1 (same sign of a on both site types: g^{−T} = 1 − aE₄₁), checked with symbolic a (cell 9). For a = 0 this is
the companion Bethe vacuum ⟨0_{w=z}|, NOT the SoV vacuum; the SoV vacuum is (λ ⊗ 1)^{⊗2}R⁻¹ = e₁⊗e₄⊗e₁⊗e₄ (checked). The overlaps
⟨W_a|Ψ_Λ⟩ = ⟨W_aR⁻¹|Ψ_c⟩ (|Ψ_c⟩ = R|Ψ_Λ⟩) were checked on two eigenstates in the notebook and on all 256 in python.

## 2. SoV expansion of ⟨W_a^L|R⁻¹ = Σ_x c_x(a)⟨x| [exact]

c_x(a) = (⟨W_a|R⁻¹ Y)_x, Y = W⁻¹ (block-diagonal in the SoV charge |x| = a₁+k₂+a₃+k₄).

| component | meaning | # nonzero words | support by charge 0..12 | site ranks | pair-cut rank (12)\|(34) | rank (13)\|(24) | rank_{Q(r₂,r₄)} of m_{a₁a₃} |
|---|---|---|---|---|---|---|---|
| a⁰ | ⟨vac| (companion Bethe vacuum ⟨0_z|) | 256 | 1,4,10,20,31,40,44,40,31,20,10,4,1 (full) | 4,4,4,4 | 11 | 4 | 2 |
| a¹ | ⟨vac|ρ(E₁₄) | 256 | full | 4,4,4,4 | 13 | 12 | 2 |
| a² | ⟨vac|ρ(E₁₄)²/2 | 256 | full | 4,4,4,4 | 16 | 16 | 4 |
| a³ | ⟨vac|ρ(E₁₄)³/6 | 256 | full | 4,4,4,4 | 13 | 12 | 2 |
| a⁴ | ⟨vac|ρ(E₁₄)⁴/24 = ⟨LW²| | 256 | full | 4,4,4,4 | 11 | 4 | 2 |
| generic a = 3/7 | full bra | 256 | full | 4,4,4,4 | 16 | 16 | 4 |

- **Degree and support.** Every c_x(a) is a degree-4 polynomial containing all five powers a⁰..a⁴ (all 256 words), so the bra is
  maximally dense in the SoV word basis, and also in every single component. (Physical frame: sparse; SoV frame: dense.)
- **z-dependence.** The words are z-independent, the bra is not: c_x(a; κz) = κ^{6−|x|} c_x(a; z) [exact, symbolic a, κ = 2/3].
- **Relations between components.** The charge grading does not separate them (all components have the full charge support and the same
  κ-scaling). The ranks are palindromic n ↔ 4−n (W_a = −a²⟨LW|ρ(e^{E₄₁/a}) per pair) but there is no simple map between the
  coefficient tensors (ratios c^{(4−n)}_{map(x)}/c^{(n)}_x for map = x→3−x, identity, A- or F-sites only: 201–256 distinct values).
  **The n = 0 and n = 4 components are exact eigenbras of t^{G_c(z)}** [exact]: n = 0 is (v(z₁)⊗ũ₄)², eigenvalue the vacuum
  τ_a(u) = Σ_I ∏_k z_{i_k}V_{i_k}(u−(k−1)h); n = 4 is (v(z₄)⊗ũ₁)², eigenvalue the same formula with (z₁,z₂,z₃,z₄) → (z₄,z₂,z₃,z₁)
  (checked t₁ at two u, proportionality exact and eigenvalues 0.0103/−0.4210 vs 0.0653/−0.3120). n = 1, 2, 3 are not eigenbras.
  So the n = 0, 4 components are SoV expansions of the vacuum and of the "lowest-weight" eigenbra and factorise on the dual side
  (§4); n = 1, 2, 3 do not.
- **Rotations by words?** Not asked here; see §6.

## 3. ⟨W_a^L|Ψ⟩/⟨0|Ψ⟩ = Σ_x c_x(a)τ_x(Ψ): exactness, sectors, structure

- **Sectors [exact + hp].** ⟨vac|ρ(E₁₄)ⁿ/n!|Ψ⟩ ≠ 0 only if Ψ has weight (2−n,0,0,n−2): sector sizes 1, 12, 28, 12, 1 (n = 0..4), 54 of the 256
  eigenstates; ⟨W_a|Ψ⟩ = aⁿ⟨W₁|Ψ⟩. The SoV sum Σ_x c^{(n)}_x τ_x(Ψ)⟨0|Ψ⟩ equals the direct overlap ⟨W^{(n)}R⁻¹|Ψ_c⟩ on all 256 states
  [hp, 50 digits: max error 5e-41, 5e-42, 5e-41, 4e-43, 5e-44 for n = 0..4], including the vanishing (≤ 4e-42) on the other 202 states.
  Since the 256 τ-monomials are independent functions on the 256 states, c_x(a) is *the* unique representation; any candidate
  formula that is valid on all 256 states is an identity of coefficient tensors.
- **A–A matrix.** ⟨W^{(n)}|Ψ⟩/⟨0|Ψ⟩ = Σ_{a₁,a₃} τ_{a₁}(θ₁−h) m^{(n)}_{a₁a₃}(r₂,r₄) τ_{a₃}(θ₃−h), r_s = τ₁(θ_s) (s = 2, 4),
  m_{a₁a₃} = Σ_{k₂k₄} c_{a₁k₂a₃k₄} r₂^{k₂}r₄^{k₄}, entries of bidegree (3,3). With τ_a(θ−h) = D_a/D_0 (dictionary, verified in §4 on vac(w)):
  rank 2 (n = 0,1,3,4) means **a sum of 2 products of bordered 4×4 Casoratians** B^{(1)}_k B^{(3)}_k/(Q₁₂₃(θ₁−h)Q₁₂₃(θ₃−h)),
  B^{(s)}_k = det[Q_i(θ_s−jh)_{i=1,2,3; j=0..3}; border row w^{(k)}_j(r₂,r₄)] (one at each A site, border rows polynomial in the F-site data r₂, r₄) —
  equivalently one 8×8 bordered determinant with two border rows (rows: 3 Casoratian rows for each A site, 2 border rows; the
  Laplace expansion over the border columns gives m_{a₁a₃} = w_{a₁}w′_{a₃} − w′_{a₁}w_{a₃}). Rank 4 (n = 2 and generic a) is no
  compression at all (any 4×4 matrix is a sum of 4 products): the singlet sector (n = 2, weight 0, 28 states) is the structurally hardest.
- **Finer structure [exact, sympy].** n = 0, 4: all entries of m share a common factor of bidegree (2,2) (a product of two bilinear
  forms in (r₂,r₄), e.g. 275968000 r₂r₄ − 145299000 r₂ + 48925184 r₄ − 44126237 and 3532390400 r₂r₄ − 3487176000 r₂ + 1174204416 r₄ − 1985680665 at the standard point), the 2×2 minors are products of such bilinear forms, and the column
  space of m is cut out by two **bilinear** equations (left kernel of bidegree (1,1)); this is the vacuum structure. n = 1, 3: no common factor of the
  entries; the 2×2 minors are products of two *irreducible* (3,3) polynomials, the column space contains a single polynomial vector of bidegree ≤ (1,1)
  (4 of bidegree ≤ (2,2), 12 of ≤ (3,3)), left kernel of bidegree (3,3): no sparser or lower-degree border structure.
- **Candidates tested on the 256 eigenstates and excluded.**
  1. Product over pairs of bordered Casoratians (pair 1 rotated × pair 2 rotated, each one-pair overlap with the other pair at the SoV
     vacuum): excluded exactly (pair-cut rank of c_x would have to be 1; it is 11–16) and numerically: G/(P₁P₂) over the 54 states with
     nonzero overlap ranges over [1.0e-6, 10.0] [hp].
  2. Linear combination αP₁ + βP₂ + γ of one-pair overlaps: LSQ residual 79% [hp].
  3. Single 8×8 bordered determinant (2 border rows): valid for n = 0,1,3,4 (rank 2) and **not** for n = 2 nor for the full bra at generic a (rank 4).
  4. Sum of a few pair-factorised terms: ≥ 11 (pair-cut rank).
  5. GLMR/2603 type Q-bracket determinant: needs a measure; the mixed 4/4̄ measure is the non-diagonal M (1456 nonzeros, triangular): not available.
  Outcome: **no pair-factorised or elegant determinant form**; the exact structure is the 256-word expansion with A–A rank 2 (n ≠ 2) as above.
- **The a¹ component alone** (⟨vac|ρ(E₁₄)|Ψ⟩, 12 states): no closed Q-form found. Exact identity (first-order perturbation theory
  in the twist, from [ρ(E₁₄), t^Λ(u)] = (z₄−z₁)T₁₄(u), not re-verified numerically): (τ_Ψ(u) − τ_vac(u))⟨vac|ρ(E₁₄)|Ψ⟩ = (z₄−z₁)⟨vac|T₁₄(u)|Ψ⟩ for all u,
  which relates it to a matrix element of the monodromy and, at special u, to few product-state overlaps, but not to Q-values.

## 4. Contrast: the twisted vacuum ⟨0_w| — an explicit finite Q-sum

- **Twist independence [exact].** Bra words ⟨x| and ket words |y⟩ built with t^{G_c(w)} equal those built with t^{G_c(z)} (256×256, exact), hence
  Gram and M = Gram⁻¹ are twist-independent. ⟨0_w| = (v(w₁)⊗ũ_w)^{⊗2} (ũ_w = row 4 of V(w)⁻¹) is an exact eigenbra of t_a^{G_c(w)}(u), a = 1,2,3, with vacuum
  eigenvalue τ_a(u) = Σ_{|I|=a}∏_{k=1}^{a} w_{i_k}V_{i_k}(u−(k−1)h) (V₁ = ∏_A(u−θ+h)∏_F(u−θ), V₂ = V₃ = ∏(u−θ), V₄ = ∏_A(u−θ)∏_F(u−θ−h)).
- **Factorisation [exact].** ⟨0_w|y⟩ = ⟨0_w|0̄⟩ ∏_s f_s(y_s), f_A(m) = τ₁(θ_A)^m, f_F(n) = τ_n(θ_F + nh), for all 256 y (relative error exactly 0);
  tensor ranks all 1. ⟨0_w|0̄⟩ = 651015225/2427364.
- **Mirror dictionary on the Q's of vac(w) [exact].** Vacuum Q-functions Q_i = w_i^{u/h}Φ(u)q_i(u), Φ(u) = V₄(u)Φ(u−h), q_i the unique minimal polynomial solution of
  Σ_a(−1)^a w_i^{−a}τ_a(u)q_i(u−ah)∏_{b=a}^{3}V₄(u−bh) = 0 (degrees 4,2,2,0); Casoratian parts q_I = det[w_i^{|I|−1−s}q_i(u−sh)]/F_{|I|} (divisions exact),
  F₂ = ∏_F(u−θ_F), F₃ = ∏_F(u−θ_F)(u−θ_F−h), F₄ = ∏_F(u−θ_F)(u−θ_F−h)(u−θ_F−2h)∏_A(u−θ_A+h). Dual Q's w_I(u) = Z_I^{(u−θ)/h}q_I(u), Z_I = ∏_{i∈I}w_i, triples I. Then
  - A site s (m-th power): τ̃ = ( det G · N₀^{A_s} · w_I(θ_s−2h)/w_I(θ_s−h) )^m, any triple I (independent of I, exact);
  - F site s (level n = 1,2,3): τ̃ = (det G)ⁿ N₀^{F_s,n} · det[w_{I_l}(θ_s+(j−1)h)]_{l=1,2,3; j∈{0..3}∖{n}} / det[…]_{j=1,2,3}, I_l = (123),(124),(134), det G = ∏w_i;
  - constants (twist- and state-independent exact rationals, N₀ = const/(det G)ⁿ; checked identical from vac(z), vac(w), vac(w₂); note N₀^{A₁} = 1/c₁, c₁ = 603680/25947 of the earlier localisation calibration):
    N₀^{A₁} = 25947/603680, N₀^{A₃} = 391645/298144;
    N₀^{F₂,n=1,2,3} = −459/19250, −26625213/5929000000, −1498387112001/130438000000000;
    N₀^{F₄,n=1,2,3} = −3589/35840, −15407577/7064780800, 229218523029/2785219182592000.
  With these, the Q-form of τ̃_y(vac_w) agrees with the eigenvalue product for all 256 y at twists w, w₂, z (error exactly 0) [exact]. The bra-side dictionary (F: τ₁(θ_F) = Q̃_j(θ_F)/Q̃_j(θ_F−h), any j; A: τ_a(θ_A−h) = D_a/D_0, any triple) also holds exactly on vac(w) and vac(z).
- **Finite Q-sum.** Completeness 1 = Σ|y⟩M_{yx}⟨x| gives ⟨0_w| = Σ_x(Σ_y⟨0_w|y⟩M_{yx})⟨x| [exact; the expansion is dense, 256 words] and

    **W_w(Ψ) := ⟨0_w|Ψ_z⟩/(⟨0_w|0̄⟩⟨0|Ψ_z⟩) = Σ_{x,y} τ̃_y(vac_w) M_{yx} τ_x(Ψ),**

  τ̃_y(vac_w) = ∏_s (mirror dictionary above on the Q's of vac_w), τ_x(Ψ) = τ_{a₁}(θ₁−h)τ₁(θ₂)^{k₂}τ_{a₃}(θ₃−h)τ₁(θ₄)^{k₄} (= D_a/D_0 Casoratian ratios of Ψ's Q's),
  M = the exact measure (1456 nonzero entries, charge-graded: M_{yx} ≠ 0 only if |x|+|y| = 12). Verified [hp] on all 256 eigenstates Ψ_z: max_n |formula/direct − 1| < 1e-44 (Mathematica, 60 digits; python 1.6e-39); min_n|W_w| = 2.9e-8 (all 256 nonzero).
  Biorthogonality as a finite Q-sum: for w = z, W_z(Ψ) = 0 on all 255 non-vacuum states (< 1e-45, relative to |W_z(vac)|), W_z(vac) = 1.3477e-4.
  No determinant form was attempted (needs the mixed-chain measure).
- **⟨W_a|y⟩ does NOT factorise [exact].** Tensor ranks of y ↦ ⟨W^{(n)}|y⟩ (flatten each site / pair cut): n = 0,4: 1,1,1,1 / 1 (eigenbras: factorised, as the twisted vacua);
  n = 1,3: 4,4,4,4 / 6; n = 2: 4,4,4,4 / 15; generic a = 3/7: site ranks 4,3,4,3, pair cut 9. Same for the bra-side: c_x has pair-cut ranks 11–16 (above).
  This is the quantitative form of "⟨W_a| shares the Q's of ⟨vac| but is not a twisted vacuum" (`w-a-reference-status.md` §3).

## 5. Global Borel rotations of the SoV vacuum ⟨0| (companion frame) [exact]

g = 1 + bE_ij on all four sites (A: v → vg, F: u → ug^{−T}); the diagonal-frame generator transported to the companion frame is g_c = V^{−T}gV^T (A), V g^{−T}V⁻¹ (F).

| generator | words at b⁰, b¹, b², b³, b⁴ (diag-frame generator transported) | words, generator acting directly in the companion frame |
|---|---|---|
| E₁₂, E₃₄, E₁₄, E₁₃, E₂₄, E₂₃ (all six) | 1, 30, 140, 236, 256 (union 256) — dense | E₁₂: 1,4,8,0,0 (13 words); E₃₄: 1,2,1 (4); E₂₄: 1,3,1 (5); E₁₃: 1,9,23 (33); E₁₄: 1,16,36,16,1 (70); E₂₃: trivial (1 word) |

So none of the physical-frame Borel rotations is few-word in the SoV basis (they are dense for the same reason ⟨W_a| is); the few-word global rotations are the E₁₂/E₃₄/E₂₄-type steps of the *companion* gl(4) (cf. README, "global rotations in the word span").

## 6. What failed / limits / open

- The shared kernel was wiped once during the session (python fallback used meanwhile); everything was re-run afterwards in the notebook.
- No Q-form for ⟨W_a^L|Ψ⟩ was found; the bra is dense in the SoV basis and pair-cut rank 11–16; the mixed 4/4̄ measure/right basis is still the missing ingredient for any 2603-type determinant.
- The minimal polynomial-border decomposition for n = 1, 3 is not free (column space has generators of bidegree (1,1) and larger); I did not construct explicit U, V with constant core.
- Only L = 2 pairs, generic θ, z (not the ABJM point); W_a with z₂ = z₃ is singular in the companion frame (README caveat).
- Not done: Q's of non-vacuum eigenstates were not recomputed in this session (the Q-dictionary on Ψ_z relies on `sov_l4.py` main(), 1e-42 on all 256 states); W_w sums use the eigenvalue form of τ_x(Ψ).

## 7. Files

`Wa-sov-expansion-L4.wb` (29 cells, all run; §6 task 4 = cells 12–21, §7 task 3 = 22–24, §8 task 5 = 25–26, 8b = 27–28, summary = 29), `checks/session-2026-10-13/wasov/`:
`wa_sov_core.py` (exact build: W, Y = W⁻¹, kets, Gram, M), `run2.py` (components, ranks, ⟨0_w|y⟩ factorisation, z-independence; `run2.log`),
`run3.py`, `run10.py`, `run11.py` (symbolic structure of m(r₂,r₄)), `run4.py` (§5, `run4.log`), `run5.py` (vacuum Q's, dual constants, mirror dictionary, `run5.log`),
`run6.py` (W_w on 256 states, `run6.log`), `run7.py` (sectors, component sums, eigenbra check, `run7.log`), `run8.py` (component relations), `run9.py` (candidate tests; writes
`Wa_sov_coefficients.json` = exact c^{(n)}_x, n = 0..4, 256 words each).

## Follow-up: W_a overlap = twisted-vacuum overlap

Notebook cells 30–48 of `Wa-sov-expansion-L4.wb` (§10.1–10.8), all run, saved with outputs; python cross-checks `checks/session-2026-10-13/wasov/run12.py`, `run13.py`.
**This supersedes the earlier statements "no closed Q-form" (§3, §6): the dense 256-word expansion is the unique representation on all 256 states, but on the 54 states of weight (μ,0,0,−μ) a 9-word F-only representation exists (below).**

**Key fact (coordinator) [exact + hp].** For w = (w₁,z₂,z₃,w₄) the diagonal-frame image of ⟨0_w| per pair is (v₁..v₄)⊗(u₁,0,0,u₄), v_j = λ_j∏_{m≠j}(w₁−z_m), u_j = ∏_{m≤3}(z_j−w_m)/∏_{m≤3}(w₄−w_m)
(checked exactly: ⟨0_w|_c R = (v⊗u)^{⊗2} with my companion normalisation). On Ψ of weight (μ,0,0,−μ): ⟨0_w|Ψ⟩ = (v₁u₄)^L⟨W_a^L|Ψ⟩, a = v₄/v₁ = −u₁/u₄ = −(z₁−z₂)(z₁−z₃)(w₁−z₁)/((z₄−z₂)(z₄−z₃)(w₁−z₄)).

**(1) Check [hp 40 digits, coordinator code in cell 31, sectors extended].** Difference |⟨0_w|Ψ⟩ − (v₁u₄)²⟨W_a|Ψ⟩| = 0 (to 40 digits) in the five sectors (2,0,0,−2) [1 state], (1,0,0,−1) [12], (0,0,0,0) [28], (−1,0,0,1) [12], (−2,0,0,2) [1] (overlaps up to 0.150, 0.073, 0.0325, 0.0114, 0.0034); both vanish in (1,1,−1,−1) [4 states]; in (0,1,0,−1) [12 states] only ⟨0_w| is nonzero (0.125), as expected. The Möbius map (cell 32, symbolic): a = v₄/v₁ = −u₁/u₄ = the formula above (difference 0), u₂ = u₃ = 0, v·u = 0, a(w₁=z₁) = 0, a(w₁=z₄) = ∞ (pole), a(w₁→∞) = κ₀ = −(z₁−z₂)(z₁−z₃)/((z₄−z₂)(z₄−z₃)) (= −245/76 at the standard point).

**(2) Full Q-chain at generic θ, z [hp 60 digits, three values of w₁ = 7/5, 11/3, −5/4 with w₄ = 5/9, 2/7, 13/9].** ⟨W_a^L|Ψ⟩ = (v₁u₄)⁻²⟨0_w|0̄⟩⟨0|Ψ⟩𝕎_w(Ψ), 𝕎_w = Σ_{x,y}τ̃_y(vac_w)M_{yx}τ_x(Ψ) with τ̃_y from the dual vacuum Q's of twist w = (w₁,z₂,z₃,w₄)
(mirror dictionary, constants N₀ of §4): relative error < 1e-46 on the 54 states for each w₁ (cell 34); diagonal image identity ⟨0_w|R == (v⊗u)^{⊗2} exact. On the other 202 states the RHS is nonzero on 72 states (those whose weight lets an A-site e₂/e₃ be compensated, e.g. (0,1,0,−1)), while ⟨W_a|Ψ⟩ = 0 there: the identity is a statement about the (μ,0,0,−μ) sectors only.

**(3) Simplification at w₂ = z₂, w₃ = z₃ [exact + hp].**
- *τ̃_y(vac_w) is independent of w₄ for all 256 y* (symbolic, cell 40) and depends on w₂,w₃ only through s = w₂+w₃, p = w₂w₃. Reason: at the ket points V₁ = 0 at θ_F, V₂ = V₃ = V₄ = 0 at θ_A, V₄(θ_F+h) = 0, V₂ = V₃. A site: τ̃ = (w₁V₁(θ_A))^m (monomial); F site: affine in w₁ (n = 3: monomial in w₁). So τ̃_y is a polynomial in w₁ of degree m₁+m₃+[n₂>0]+[n₄>0] ≤ 8.
- Hence 𝕎_w(Ψ) is a **polynomial of degree 8 in w₁** and on the (μ,0,0,−μ) sectors (n = 2−μ) it is exactly
  **𝕎_w(Ψ) = [λ₁(z₄−z₂)(z₄−z₃)(w₁−z₂)(w₁−z₃)(w₁−z₄)²]² a(w₁)ⁿ Γ_Ψ,  Γ_Ψ = ⟨W₁|Ψ⟩/⟨0|Ψ⟩**, i.e. ∝ (w₁−z₂)²(w₁−z₃)²(w₁−z₄)^{4−n}(w₁−z₁)ⁿ
  (w₄ cancels against ⟨0_w|0̄⟩ = 1/∏_{m≤3}(w₄−w_m)²) [hp, error < 1e-46, three (w₁,w₄), cell 48]. So all w₁-dependence is this explicit factor: one number Γ_Ψ per state, and the Q-sum at any w₁ gives it.
- *Which (x,y) terms survive.* The coefficient of w₁^j in Σ_yτ̃_yM_{yx} (a vector over x) has 9, 32, 64, 96, 110, 96, 64, 32, 9 nonzero words for j = 0..8 (palindromic). The top power w₁⁸ comes only from the 9 ket words y = (3,n₂,3,n₄), n₂,n₄ ∈ {1,2,3}, and is supported on exactly the 9 **F-only** words (0,k₂,0,k₄), k₂,k₄ ≤ 2, and **equals the SoV expansion of e₁⊗ũ⊗e₁⊗ũ, ũ = (0,1,−s,p) (coefficients of (x−z₂)(x−z₃)), with ratio exactly 1** (companion frame; w₁→∞: v(w₁) → w₁³e₁ and u_w → ũ/((w₄−z₂)(w₄−z₃))). By the identity the lower coefficients carry no independent information on the sectors.
- **Compact formula [hp 1e-46 on 54 states; coefficients exact; general-θ form verified exactly at two (h,θ) points].** For Ψ of weight (2−n,0,0,n−2):

  ⟨W_a^L|Ψ⟩ / ⟨0|Ψ⟩ = (a/κ₀)ⁿ · Σ_{b,b′∈{2,3,4}} ũ_b ũ_{b′} P_{bb′}(ρ₂,ρ₄) / (λ₁²((z₄−z₂)(z₄−z₃))²),  ũ = (0,1,−(z₂+z₃),z₂z₃)  (b, b′ over components 2,3,4 of ũ),

  ρ₂ = τ₁(θ₂)/(hθ₁₂θ₂₃(θ₂₄−h)), ρ₄ = τ₁(θ₄)/(hθ₁₄θ₂₄θ₃₄), τ₁(θ_F) = Q_j(θ_F)/Q_j(θ_F−h) (single-Q ratio, any j; Γ gauge of sov_l4), ν = θ₂₄, κ₀ as above, λ₁ = 1/∏_{m≠1}(z₁−z_m), and
  P₄₄ = 1, P₄₃ = ρ₄ − (h/ν)ρ₂, P₃₄ = ρ₂, P₃₃ = (ν/(ν+h))ρ₂ρ₄, P₄₂ = ρ₄² − (h/ν)(1+h/(ν+h))ρ₂ρ₄ − (h/ν)((ν−h)/ν)ρ₂², P₂₄ = ((ν−h)/ν)ρ₂² + (h/(ν+h))ρ₂ρ₄,
  P₃₂ = (ν/(ν+h))ρ₂ρ₄² − (h/(ν+h))ρ₂²ρ₄, P₂₃ = (ν/(ν+h))ρ₂²ρ₄, P₂₂ = (ν/(ν+h))²ρ₂²ρ₄² (P_{bb′} = coefficient for F-site vector e_b at site 2, e_{b′} at site 4; the e₃/e₄ block is the reports §6 formula).
  So **⟨W_a^L|Ψ⟩ at L = 2 is a polynomial of bidegree (2,2) in the two F-site single-Q ratios τ₁(θ₂), τ₁(θ₄)** (times ⟨0|Ψ⟩ and explicit factors); only the 9 words a₁ = a₃ = 0, k ≤ 2 are needed instead of 256.
  Verified: direct vs formula at generic θ and at the ABJM pair-by-pair θ = (2/5,−3/5,−3/8,−11/8) with generic z: error < 1e-51 on all 54 states (cell 42); the coefficients read off at the standard point are reproduced exactly at (h,θ) = (2/3,(3/11,−5/7,9/13,−13/9)) (run12.py). The formula involves only θ-differences with ν = θ_F1−θ_F2 ≠ 0, ±h and θ₁₂ = θ₃₄ = 2h: **regular at the ABJM point** (pair-by-pair θ).
- *Negative/limiting remarks.* (i) The 256-word expansion c_x(a) of §2–3 is still the unique representation valid on all 256 states (it encodes the vanishing off-sector); the 9-word form is valid only on the 54 sector states, which is why the dense expansion looked structureless. (ii) Tested for L = 2 only; for L pairs the same mechanism gives an F-only bra with L F-sites (words k_s ≤ 2) — not tested. (iii) ⟨0|Ψ⟩ is still needed: the formula expresses ⟨W_a|Ψ⟩ through Q-ratios relative to ⟨0|Ψ⟩; its absolute normalisation (SoV norm 𝒩 = 1/Στ̃Mτ at the ABJM point, regularised bases) was not rederived. (iv) The vacuum Q's/dual Q's themselves were not found to simplify at w₂=z₂, w₃=z₃ beyond the τ̃ statement above.

**(4) ABJM point [hp].**
- *Identity.* The twisted-vacuum identity is algebraic in the diagonal frame; verified on the pair-by-pair ε-chain (z = 1+εφ): (v₁u₄)⁻²⟨0_w|Ψ_ε⟩ = ⟨W_a|Ψ_ε⟩ to 1e-83 on all 28 states of the zero-weight sector at ε = 1/100 and 1/400 (cell 46), with w₁ = (a z₄ − κ₀ z₁)/(a − κ₀) = 1+εω₁, ω₁ = (aφ₄−κ₀φ₁)/(a−κ₀), κ₀ → −(φ₁−φ₂)(φ₁−φ₃)/((φ₄−φ₂)(φ₄−φ₃)) (ω₁ = −0.7351575456 for φ = (0.9,0.35,0.6,−1.1), c = 0.3, a = 0.3). I did **not** use the regularised bra/ket bases: the closed formula needs only the F-only words ⟨0|t₁(θ₂)^{k₂}t₁(θ₄)^{k₄}, which are regular at the ABJM point, so no regularisation was required for ⟨W_a|Ψ⟩.
- *Even bridge M = 1 at L = 2 through Q-data (cell 44).* z = 1+εφ, φ = (9/10,7/20,3/5,−11/10) (φ₂ ≠ φ₃; K_κ with y = (z₂+z₃)/2), θ = (θ_F1+2h,θ_F1,θ_F2+2h,θ_F2), θ_F = ε(0,0.37), c = 0.3, κ = 1−cε, a = (φ₁−c)/(φ₁−φ₄), unit-norm Ψ_ε, exact rational operators, 90-digit eigenvectors, ε = 1/100, 1/200, 1/400, 1/800:
  C_Q(ε) = K_κ(τ_Ψε) · [⟨0|Ψ_ε⟩(a/κ₀)²ũ·P(ρ₂,ρ₄)·ũ/(λ₁²((z₄−z₂)(z₄−z₃))²)] / (1−κz₁)²,
  K_κ = s²+κs[X₁/(h²λ₁)+X₂/(2h⁴yλ₂)+X₃/(6h⁵y²λ₃)]+κ²X₃X₁/(24h⁷y²λ₁λ₃) from τ₁(θ₂), τ₂(u₀), reduced τ₃(u₀) (pair-2 localisation constants λ_k). C_Q = C_direct = K⟨W_a|Ψ⟩/N to |w−w_Q|/|w| < 1e-66 for all 28 states and all ε. The two singlets:
  |C|(ε) = 0.1579566, 0.1580350, 0.1580744, 0.1580941 → cubic fit **0.15811388301 vs 1/√40 = 0.15811388301 (diff 2e-12)**; 0.2039019, 0.2040176, 0.2040720, 0.2040983 → **0.20412414524 vs 1/√24 = 0.20412414523 (diff 9e-12)**. (Double-precision python `run13.py` reproduces the formula at ε ≥ 0.1 only; the cancellation in ⟨0|Ψ⟩ ~ Σλ_jψ_j with λ ~ ε⁻³ needs the hp of the notebook.)
- *What is Q and what is not.* In C_Q: K_κ (Q-forms of X₁, X₂, X₃ from lpairs Task 4), ρ₂, ρ₄ (single-Q ratios) and explicit θ, z factors are Q-data. Not Q-expressed here: ⟨0|Ψ_ε⟩/‖Ψ_ε‖ (SoV normalisation; at the ABJM point it needs the regularised Gram/measure of reports §9d), taken from the eigenvector. So the even bridge is "Q-complete up to the norm". The Q-sum 𝕎_w(Ψ) with the twisted vacuum at w₁ = 1+εω₁ gives the same ⟨W_a|Ψ⟩ (cell 46) but needs the companion frame (z₂ ≠ z₃) and the Gram matrix M, which degenerates as ε → 0; the closed 9-word formula does not.
