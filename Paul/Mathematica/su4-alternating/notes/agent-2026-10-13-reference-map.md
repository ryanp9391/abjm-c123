# Which references admit the pair κ-family? The twisted vacuum ⟨0_w| and the map of the Borel orbit (2026-10-13)

Agent report (research agent, 2026-10-13). Notebook: `reference-admissibility-Kkappa.wb` — **all sections run in the Wolfram kernel** (sections 0–5 before the
shared kernel died, everything else after its restart at 08:03). While the kernel was monopolised/dead (>10 min) I wrote an exact python **fallback**
(`checks/session-2026-10-13/refmap/`, `Fraction`/`sympy`/mod p on top of `pgbase.py`/`eng.py`); every Wolfram output of §5c–§9c was afterwards compared with it and agrees
(notebook §10). Items marked "python only" below were not ported to Wolfram (exhaustive {−1,0,1,2}⁴ scan, symbolic coefficient families for classes R5/R2, `K^F` on the chain).
Evidence tags: [exact] rational / symbolic; [mod p] exact mod p = 268435399 (false-positive probability ~1e-9 per component); [exact→float] exact
rational least squares printed as floats.
Conventions as in README / `w-a-reference-status.md`: A F A F …, bras row vectors, A: v → v g, F: u → u g^{−T}, ρ(E_ij) = E_ij ⊗ 1 − 1 ⊗ E_ji;
ABJM point θ = (θ_F+2h, θ_F), h = 1/2; Λ = diag(z); pair-level words X₁ = t₁(θ_F), X₂ = t₂(θ_F+h), X₃ = t̂₃(θ_F+h); one-pair chain θ = (1, 0).

## 0. Result in five lines

1. **Twisted vacuum ⟨0_w|.** At the ABJM-point twist z₂ = z₃ its diagonal-frame image degenerates: v → (0, 1, −1, 0) (A site **frozen** on the gl(2)-odd
   direction, projectively independent of w), u = u_w generic. On this limit reference a common pair rotation exists [exact, L = 2, 3 mod p], but it is an **F-site
   rotation only** and is *not* the five-word K_κ of ⟨W_a|: the A site is an eigenvector of Λ, so "rotation" means u → u(1−κΛ)^{-1}, v → v·const. The operator is
   K^F_κ = [1 − κe₁′ + κ²e₂′] + 2κ X₁ + 8κ² X₂ (pair level; e₁′ = z₁+y+z₄, e₂′ = yz₁+yz₄+z₁z₄; X₃, X₃X₁ not needed), the one-site analogue of 2603's word 1 − κt₊ + κ²t₀.
2. **At z₂ ≠ z₃ no κ-family exists on ⟨0_w|** (generic w₁ ∉ {z}): span ∩ (span of the 13 κ-components) has dimension 2 of 13; the K-family exists iff w₁ ∈ {z₁,…,z₄}
   (then v = e_a frozen). Near z₂ = z₃ (z₂−z₃ = εΔφ) the best-fit residual of T_κ is **linear in Δφ and charged** (§6), so it does not heal in the untwisting scaling.
3. **Admissible set (closed form, exact).** At z₂ = z₃: K-family exists ⇔ **[v₁u₁ = 0, v₄u₄ = 0, v₂u₂ + v₃u₃ = 0]** (block-BPS: v|_C·u|_C = 0 in every spectral block of Λ)
   **or [v₂ = v₃ = u₂ = u₃ = 0]** (both sites on indices {1,4}: the torus images of W_a). At generic z only the block condition (v_j u_j = 0 for every j) survives; W_a fails.
   Checked on all 8946 references with v,u ∈ {−1,0,1,2}⁴, v·u = 0 (0 mismatches), 196 structured random samples, and by symbolic `Reduce` over all Levi classes.
4. **Why only W_a is special / no-go.** The cross-block pairings that make the two-block polarisation sensitive to s > 0 states are p₂₃ = Σ π_C/(c−φ_C), p₃₂ = Σ π_C(c−φ_C),
   π_C = v|_C·u|_C. Every admissible reference except the W-type has all π_C = 0 ⇒ **p₂₃ ≡ p₃₂ ≡ 0**; W_a is admissible exactly because π₁ = −π₄ ≠ 0 is allowed at z₂ = z₃.
   So (cross-pairing ≠ 0) and (K-admissible) are compatible only for the W-type family.
5. **Untwisting.** κ = 1 − cε, z = 1+εφ: 1 − κz_j = ε(c − φ_j): **D = diag(c − φ_j), not c + φ_j** [exact, ε = 10⁻²…10⁻⁵]. Generic φ: T_κ → (v̂D)⊗(ûD⁻¹) with the stated v̂, û (deviation O(ε)); but
   at φ₂ = φ₃ the same formula collapses to (0,1,−1,0) ⊗ (û_j/(c−φ_j)): A frozen, p₂₃ = p₃₂ = mv = a₂ = a₃ = 0.

## 1. Set-up actually used (pair level, with the L = 2, 3 verification)

By the column-support lemma every localising word on (v⊗u)^{⊗L} equals λ_w × (same word on the 1-pair chain) ⊗ (v⊗u)^{⊗(L−1)}; re-verified here at L = 2 and **L = 3 pairs (6 sites, 4096-dim) [mod p]**
for every reference of §3/§4 that passes (and for the negative controls): localisation of all five words, and
Σ_w (c_w/λ_w) ψ·word_w = T_κ/(1−κy)² ⊗ (v⊗u)^{⊗(L−1)} for κ ∈ {3/10, −7/5, 1/z₄}, c_w the pair-level solution (free ones = 0)
(`l3check.py`, `l3check_out.txt`; Wolfram version: notebook §7, same result).
Target: T_κ = [v(1−κΛ)] ⊗ [u adj(1−κΛ)] (polynomial normalisation; = ∏_m(1−κz_m) × the group-rotated vector).
Calibration: for W_a the pair-level coefficients are exactly K_κ's: c = (s², κs/h², κs/(2h⁴y), κs/(6h⁵y²), κ²/(24h⁷y²)), s = 1−κ(z₁+z₄) [exact, sympy, symbolic in z₁, y, z₄, κ]; at κ = 3/10, z = (3/2,5/4,5/4,4/7) these are
(2809/19600, 159/350, 636/875, 1696/4375, 192/625), identical to the Wolfram calibration cell.

**Charged-sector structure [exact, symbolic z].** The words commute with the Cartan torus; on the 12 charged states e_a⊗e_b (a≠b) with x₁ = X₁/2h², x₂ = −X₂/2h⁴, x₃ = X₃/12h⁵ the eigenvalues are
x₁ → z_b, x₂ → z_b(σ − z_a − z_b), x₃ → det/z_a, x₃x₁ → det z_b/z_a (σ = Σz, det = Πz). So on a reference with charged support S = {(a,b): a≠b, v_a u_b ≠ 0} the five words act through the |S|×5 evaluation matrix; a κ-family needs
each κ-component (a unit-like vector on a class pair) to lie in the span, i.e. rank |S| (and a compatible zero-weight sector). At generic z, S has 12 elements for a generic reference: impossible.

## 2. Calibration and the four references at z₂ = z₃ (scan items (1)–(4)) — Wolfram cell 11 and python `scan.py`, identical

z = (3/2, 5/4, 5/4, 4/7), κ ∈ {3/10, −7/5, 1/z₄}; "K" = T_κ ∈ span for all three κ.

| reference (v ; u) | rank span | K-family |
|---|---|---|
| (1) W_a, a = 1/3, −1: (1,0,0,a) ; (−a,0,0,1) | 5 | **yes** (coefficients = K_κ's) |
| (1) bare vacuum e₁ ; e₄ | 1 | yes |
| (2) random BPS (1,2,−1,3) ; (5/3,−1/6,4/3,0) | 5 → 6 | no (κ = 1/z₄ only: rank 5) |
| (2) random BPS (1,−2/3,5/7,3/2) ; (−19/21,1/4,−3/5,1) | 5 → 6 | no |
| (3) ⟨0\|_comp image, **formal** (generic z₃ = 2/3 formula used on the z₂=z₃ chain) | 5 → 6 | no |
| (3) ⟨0\|_comp image, true z₂=z₃ limit: (0,1,−1,0) ; (1,1,1,1) | 3 | **yes** (not with K_κ's coefficients, see §4) |
| (4) ⟨0_w\| formal | 5 → 6 | no |
| (4) ⟨0_w\| z₂=z₃ limit, generic w: (0,1,−1,0) ; u_w = (9/224, −69/4480, −69/4480, 899/13720) | 3 | **yes** |
| (4) locus w₁ = z₁: e₁ ; (0,−23/896,−23/896,403/5488) | 2 | yes |
| (4) locus w₁ = z₄: e₄ ; (585/1568,437/6272,437/6272,0) | 2 | yes |
| (4) locus w₂ = w₃ (limit) | 3 | yes |
| (4) locus w₂ = w₃ = y (limit): (0,1,−1,0) ; (u₁,0,0,u₄) | 2 | yes |
| (4) locus w₄ = z₄ (u independent of w₄) | 3 | yes |
| (4) locus w₁ = z₁ & w₂ = w₃ = y: e₁ ; e₄ | 1 | yes |
| (4) locus w₁ = y: (0,α,β,0) ; (u₁,0,0,u₄) | 2 | yes |

"formal" rows in the Wolfram scan mix a generic-z reference with z₂=z₃ words and are only controls; the consistent generic-z test is §6.
Every "yes" among the twisted-vacuum rows has the A site frozen (v ∈ span{e₁}, {e₄}, or the Λ-eigenvector (0,p) with p·q = 0): no twisted-vacuum row has an A site with two independent directions in the sense K_κ needs.

## 3. The admissible set at z₂ = z₃ (scan item (5)) — symbolic `Reduce` over all Levi classes (Wolfram cell 13) + closed form

The centraliser of Λ is GL(1)×GL(2)×GL(1) (commutes with all t_a and with 1−κΛ), acting v → vk, u → uk^{−T}; every v is Levi-equivalent to one of seven representatives and u is kept symbolic.
Write p = (v₂,v₃), q = (u₂,u₃). Quantifier elimination (conditions on u in Levi coordinates; "True" = every u with v·u = 0):

| Levi class of v | admissible ⇔ |
|---|---|
| (1,1,0,1) (v₁,p,v₄ ≠ 0) | u₂′ = 0 and u₄ = 0 (⇒ u₁ = 0): u = (0, p^⊥, 0) |
| (1,1,0,0) (v₁,p ≠ 0, v₄ = 0) | u₂′ = 0 (⇔ p·q = 0, then u₁ = 0): u = (0, q, u₄), q ⟂ p |
| (0,1,0,1) | u₄ = 0 (⇒ p·q = 0) |
| (1,0,0,1) (p = 0) | q = 0 **or** u₄ = 0 |
| (1,0,0,0) | always (u₁ = 0 forced) |
| (0,1,0,0) | always (p·q = 0 forced) |
| (0,0,0,1) | always (u₄ = 0 forced) |

Collecting: **𝒜 = 𝒜_block ∪ 𝒜_W**,
- 𝒜_block = { v₁u₁ = 0, v₄u₄ = 0, p·q = 0 } (hence v·u = 0), i.e. v|_C·u|_C = 0 for C = {1}, {2,3}, {4}; equivalently u f(Λ) vᵀ = 0 for **every** function f;
- 𝒜_W = { p = 0 = q } (v, u supported on {1,4}, v₁u₁ = −v₄u₄ ≠ 0 allowed): the torus images of W_a.

Inside the Borel orbit (v₁u₄ ≠ 0): **𝒜 ∩ Borel = (C) ∪ (B)**, with (C) = {v₄ = 0, u₁ = 0, p·q = 0}: v = (1,p,0), u = (0,q,u₄), and (B) = {p = q = 0}: v = (1,0,0,v₄), u = (−v₄u₄,0,0,u₄)/1 (v₁u₁ + v₄u₄ = 0), i.e. W_a up to torus.
(C) contains (A) = {v = e₁, u = (0,q,u₄)} and the A-site-non-frozen references v = (1,p,0): these are the vacua rotated by g = exp(p·E_{1,(23)}) exp(c·E_{(23),4}) with p·c = 0
(the "E₁₂/E₃₄-type" good rotations of the README); they are *good* (uΛvᵀ = 0), W_a is not (uΛvᵀ = z₁v₁u₁ + z₄v₄u₄ = (z₁−z₄)v₁u₁ ≠ 0). Block-BPS ⊂ good (it kills Σ_C z_C π_C and Σ_C π_C together), the converse fails.
Closure points (v₁ = 0 or u₄ = 0) with the same conditions also pass.

Named sub-families (python and Wolfram notebook §5c, identical, exact):
- v = (1,0,0,v₄), u generic: admissible only for q = 0 (W) [or u₄ = 0, outside the orbit].
- v generic (1,p,v₄ ≠ 0), u = (−v₄,0,0,1): **fails**; gl(2)-symmetric (1,v₂,v₂,v₄) ; (u₁,u₂,u₂,1): fails unless v₄ = 0, u₁ = 0 and v₂u₂ = 0.
- only v₂,v₃ nonzero (any u with p·q = 0): **always admissible**; only u₂,u₃ nonzero: **always admissible**.
- v = (1,2,3,0), u = (0,3,−2,5): admissible (rank 4); v = (1,2,3,0), u = (−5,1,1,3) (p·q = 5): fails; v = (1,2,3,2), u = (−2,3,−2,1) (p·q = 0 but v₄u₄ ≠ 0): fails.

Exhaustive and random confirmation [exact]: 8946 pairs (v,u) ∈ {−1,0,1,2}⁸ with v·u = 0: 2690 admissible, 0 mismatches with the closed form (`smallint.py`); 39 blockwise-by-construction + 20 W-type + 137 random-zero-pattern
samples: predicted = observed in all (`scan.py`).

**Generic z (z₂ ≠ z₃)** [exact, `genericz.py`, `smallint.py`]: enumerating all zero patterns with generic values: the admissible patterns are exactly those with **disjoint supports of v and u**
(50 of 117 patterns); in the small-integer scan 1890/1890 blockwise candidates (v_ju_j = 0 ∀ j) are admissible and **no** non-blockwise reference is. W_a (and all W-type) fail at z₂ ≠ z₃.

## 4. The operators are class-dependent (K_κ is not universal)

- K_κ with the W_a coefficients acts as the rotation only on 𝒜_W (checked on W_a; on the comp-vacuum limit, on R5, R6, R2 it does **not**) — this explains the old "K_κ fails on ⟨0|_comp" (reports §7, 9d): true for these coefficients, and at generic z.
- Symbolic solutions of T_κ ∈ span (sympy, exact):
  - A-frozen class (0,1,−1,0)⊗(u₁,u₂,u₂,u₄) (= twisted vacuum limit, companion-vacuum limit, Levi image of v = e₂ with p·q = 0): c₃ = 8κ² fixed, c₁ = 1 − κe₁′ + κ²e₂′ − (3/8)c₄ y z₁z₄, c₂ = 2κ − (3/8)c₅ y z₁z₄, c₄, c₅ free; **u-independent**. Choice c₄ = c₅ = 0:
    **K^F_κ = [1 − κ(z₁+y+z₄) + κ²(yz₁+yz₄+z₁z₄)] + 2κX₁ + 8κ²X₂**; on an F index b it acts as ∏_{j∈{1,y,4}}(1−κz_j)/(1−κz_b) = 1 − κ(e₁′−z_b) + κ²(e₂′ − z_b(e₁′−z_b)). Chain level: c₁ + (2κ/λ₁)X₁ + (8κ²/λ₂)X₂.
    Verified on the chain, **L = 2, 3 [mod p]**, on the twisted-vacuum limit, the companion-vacuum limit and a generic (u₁,u₂,u₂,u₄), κ ∈ {3/10, −7/5, 7/4, 2/9} (`r6_l3.py`).
  - Class (v = e₁, u = (0,u₂,u₃,u₄)): 3 free coefficients; class (1,1,0,0)⊗(0,0,u₃,u₄): one free coefficient, other coefficients **rational in κ with a pole at κ = 1/y** (unlike the polynomial K_κ).

## 5. GL(4) covariance (item 6) [exact, `covariance.py`]

- ρ(g) t_k^Λ ρ(g)^{−1} = t_k^{G′} for k = 1,2,3 with **G′ = g^{−T} Λ gᵀ = (gΛg^{−1})ᵀ** (the other three candidates g^{−1}Λg, gΛg^{−1}, gᵀΛg^{−T} fail), hence ⟨vac|ρ(g) X[t^Λ] = ⟨vac| X[t^{G′}] ρ(g) for every word, verified exactly for the whole operator K for seven g's (W_a, torus×W_a, R5, R2, a broken R2, a generic Borel element, W_a×R2).
- For upper-triangular g, G′ is **lower** triangular with spectrum z, and ⟨vac^L| is its vacuum eigenbra (for W_a: G′ = Λ + c E₄₁, as in `w-a-reference-status.md`).
- Since T_κ = ⟨vac|ρ(g(1−κΛ)) = ⟨vac|ρ(1 − κG̃)ρ(g), G̃ = gΛg^{−1} = G′ᵀ, "common rotation on ⟨vac|ρ(g)" ⇔ ⟨vac|K[G′] ∝ ⟨vac|ρ(1−κG′ᵀ)  (verified for the W-type rows with the W_a coefficients).
- u f(Λ) vᵀ = [f(G′)]₄₁ for v = e₁g, u = e₄g^{−T} (checked: the column "[G′ⁿ]₄₁" vanishes exactly on the block-admissible rows). Hence
  **admissible ⇔ (i) [G′ⁿ]₄₁ = 0 for n = 1,2,3 (⇔ e₄ f(G′) e₁ᵀ = 0 for all polynomials f; for lower-triangular G′: G′₄₁ = 0, G′₄₂G′₂₁ + G′₄₃G′₃₁ = 0, …) or (ii) G′ = (Levi/torus conjugate of) Λ + cE₄₁, i.e. exactly one nonzero off-diagonal entry, at (4,1).**
  Rows: R5 (G′ off-diagonal at (4,2),(4,3) only) and R2 (G′ at (2,1),(3,1),(4,2),(4,3) with G′₄₂G′₂₁ + G′₄₃G′₃₁ = 0) satisfy (i); W_a, torus×W_a satisfy (ii); the broken R2 (G′²)₄₁ = −4845/784 ≠ 0, the generic Borel element and W_a×R2 fail both.

## 6. Obstruction analysis for the twisted vacuum at z₂ ≠ z₃ [exact / exact→float]

Consistent generic-z test (words and reference at the same z = (3/2,5/4,2/3,4/7), w = (7/5,9/8,3/7,5/9)), `obstruction.py`:
rank of the 5-word span 5; the 13 κ-components D, B_ab have rank 13; **dim(span ∩ span of κ-components) = 2** (so 11 of the 13 needed directions are missing); charged support |S| = 12 with the 12×5 evaluation matrix of rank 5; rank(span + T_κ) = 6 at all three κ. (W_a at the same generic z: |S| = 2, charged part fine, zero-weight leftover, rank 6.)
A "6th localised operator" cannot help: the missing space has dimension 11 and the localised algebra on pair 1 at the ABJM point is 5-dim (words of length ≤ 4, RESULTS §3).

Scaling regime z = 1+εφ, w = 1+εω, κ = 1−cε, φ = (9/10, 7/20, 7/20+Δφ, −11/10), ω = (7/10, −3/10, 1/5, 1/2), c = 4/5: relative best-fit residual of T_κ in the 5-word span
(exact least squares → float; ε = 10^{−e}; identical to the 80-digit Wolfram result where run):

| Δφ | e | tw. vac, all comps | tw. vac, **charged** comps | W_a, all | W_a, charged |
|---|---|---|---|---|---|
| 0.01 | 2/3/4 | 0.0572/0.0574/0.0574 | 0.0266/0.0270/0.0271 | 0.0113/0.0115/0.0115 | 0 |
| 0.1 | 2/3/4 | 0.342/0.343/0.343 | 0.160/0.161/0.161 | 0.112/0.114/0.114 | 0 |
| 0.5 | 2/3/4 | 0.495/0.455/0.451 | 0.337/0.309/0.307 | 0.486/0.493/0.494 | 0 |

Both residuals are **linear in Δφ = (z₂−z₃)/ε at small Δφ and independent of ε** (no healing as ε → 0). For W_a the leftover is purely in the zero-weight sector (E₂₂, E₃₃ directions; the known z₂≠z₃ effect of §4 of the 10-12 reports), for the twisted vacuum it contains **charged** components (relative ~0.16 at Δφ = 0.1).
Whether the zero-weight leftover is harmless for the physical overlaps (as observed for W_a, reports 9c) is untested for the charged leftover; there is no reason for it to be suppressed.
Quantitatively the block-BPS violation is v₁u₁/(v₂u₂), v₄u₄/(v₂u₂) = O(Δφ) because v̂₂ = −v̂₃ ~ 1/Δφ.

## 7. Lie algebra of pair-1 rotations inside the span (item d) [exact, `covariance.py` §9]

𝔥_r = {X ∈ gl₄: r ρ(X) ∈ span of the five words}; stabiliser dimension 10, orbit tangent dimension 6 for all references below.

| reference | rank span | dim 𝔥_r/stab | representatives mod stabiliser |
|---|---|---|---|
| W_a (a = 1/3) | 5 | 2 | E₁₁, E₁₄ (torus and the a-shift) |
| R2: v=(1,2,3,0), u=(0,3,−2,5) | 4 | 3 | E₁₁, E₁₂ − (3/5)E₁₄ + (1/2)E₃₃ + (1/5)E₃₄, E₂₂+E₃₃ |
| R5: v=e₁, u=(0,2,3,5) | 2 | 2 | E₁₁, E₂₂ + (3/5)E₃₄ |
| twisted-vacuum limit (0,1,−1,0) ; u_w | 3 | 3 | E₁₁, E₂₁+E₃₁ (F-site mixing u₁ → e₂+e₃), E₂₂+E₃₃ |
| R6: (0,2,3,0) ; (1,3,−2,4) | 3 | 3 | E₁₁, E₂₁ − (2/3)E₃₁, E₂₂+E₃₃ |
| random BPS (not admissible) | 5 | 1 | one mixed element |

Torus-orbit test: for **every** admissible reference the whole torus orbit {[v f(Λ)] ⊗ [u f(Λ)^{−1}]} for f = Λ, Λ^{−1}, 1−κΛ^{−1}, Λ², generic diagonal (f₂ = f₃) lies in the span (each κ-component is in the span, §5d); for the non-admissible random BPS only f = Λ^{−1} does (the single word X₃X₁ ~ ρ(Λ^{−1}), RESULTS §3).
So the 5-dim span is, on W_a, "r, ∂_t(torus), ∂_a (E₁₄) + two non-group directions"; on the A-frozen class it is the full torus orbit span plus the mixing E₂₁+E₃₁.

## 8. Untwisting [exact, `untw.py`]

z = 1+εφ, w = 1+εω, κ = 1−cε (so 1−κz_j = ε(c−φ_j) + O(ε²)); φ = (9/10,7/20,7/10,−11/10) (φ₂≠φ₃), ω = (7/10,−3/10,1/5,1/2), c = 4/5:

| ε | deviation of T̂_κ from (v̂D)⊗(ûD⁻¹), D = c − φ | D = c + φ |
|---|---|---|
| 10⁻² … 10⁻⁵ | 3.0e−2, 3.0e−3, 2.9e−4, 2.9e−5 (∝ ε) | 1.0, 0.98, 0.98, 0.98 (wrong) |

with v̂_j = ∏_{m≠j}(ω₁−φ_m)/(φ_j−φ_m), û_j = ∏_{m≤3}(φ_j−ω_m): the conjectured two-block form is confirmed with **D_j = c − φ_j** (the "c+φ" of the brief corresponds to κ = 1 + cε or φ → −φ).
At φ₂ = φ₃ (7/20) with the limit reference (0,1,−1,0)⊗u_w: deviation from (0,1,−1,0) ⊗ (û_j/(c−φ_j)) = 3.0e−2, 3.0e−3, 3.0e−4, 3.0e−5. The normalised v̂ as φ₃ → φ₂ (d = 0.1, 0.01, 0.001): (0.148,−0.524,0.839,0.002), (0.016,−0.693,0.721,0.0003), (0.0016,−0.706,0.709,0.00003) → (0,−1,1,0)/√2.
Two-block Borel invariants of the A-frozen polarisation (blocks (v,u), (vD,uD⁻¹)): a₂ = a₃ = 0, m_v = 0, p₂₃ = Σπ_C/(c−φ_C) = 0 = p₃₂; only d₂, d₃, m_u survive. With the structure of reports §8 (F polynomial in a₂,a₃,d₂,d₃,m_v,m_u,p₂₃,p₃₂ with e₈+e₉ = s) this polarisation can only see the e₈ = e₉ = 0 monomials; **it is not the generic two-block polarisation** that the 2603-type construction was meant to produce [inference from §8 of the 10-12 reports, not recomputed].
For the A-non-frozen admissible family (C) the same algebra gives a₂ = 1, d₂ = u₄, m_v = p₂(φ₁−φ₂), m_u ≠ 0 but again p₂₃ = p₃₂ = 0.

## 9. Corrections to earlier statements and caveats

- `w-a-reference-status.md` §2, last remarks: "Good references (⟨0|, the locus) do not admit K_κ" — precise statement: the **five-word K_κ of W_a** is not the rotation on them; at z₂ = z₃ the good references 𝒜_block *do* admit a (different, class-dependent) κ-operator (§4), but those references have p₂₃ ≡ p₃₂ ≡ 0.
- The diagonal-frame image of ⟨0_w| exists only for z₂ ≠ z₃ (V(z) singular at z₂ = z₃); at z₂ = z₃ it is a limit (rescaled by δ): v = (c/δ)(0,−1,1,0) + O(1), the O(1) part is lost; the exact statement is for the rescaled limit reference.
- The Wolfram "formal" scan rows (generic-z reference on the z₂ = z₃ chain) are controls, not the consistent test.
- The shared Wolfram kernel died during this session (another agent's runaway memoised cell) and was restarted by the coordinator; all my `ra*` definitions were lost and the notebook was re-run from the top (no memoisation in the new cells; every cell < 15 s).
  The python scripts were written during the outage and agree with the Wolfram outputs (rows, ranks, coefficients, residuals, L = 3 mod p, covariance, Lie algebra, generic-z patterns, symbolic K^F).
- Twisted vacuum at generic z, loci in w (Wolfram §9b and python `twloci_generic.py`, exact): only **w₁ ∈ {z₁,…,z₄}** (then v = e_a, π_j = v_ju_j = 0 for all j) admits a κ-family; generic w, and w₂,w₃ ∈ {z} (π_j ≠ 0 for two j's) do not.
- Not done: the physical test of the charged leftover at z₂−z₃ = O(ε) (overlaps with Ψ_ε), the generic-z admissible set for non-generic values beyond the {−1,0,1,2}⁴ scan, W-type and block families for other twists (only z = (3/2,5/4,5/4,4/7) and generic (3/2,5/4,2/3,4/7) used; the closed form is polynomial in z but only tested at those).

## 10. Consequences for the even-bridge program

- Route A of `w-a-reference-status.md` §4 (change the reference to ⟨0_w|, 2603's route) does **not** reproduce the 2603 structure at the ABJM point: with z₂ = z₃ the reference collapses to an A-frozen one, with z₂ ≠ z₃ the five localised words cannot rotate it (11 of 13 κ-directions missing, charged residual ∝ z₂−z₃).
- The W_a family stays the only K-admissible reference with a nonzero cross-block pairing (§0.4); every other admissible reference gives p₂₃ = p₃₂ = 0.
- If an A-frozen polarisation is acceptable (it is blind to s > 0 states, §8), the usable κ-word is the two-term K^F_κ of §4 with eigenvalue P^F_κ(τ) = 1 − κ(e₁′ − x₁) + κ²(e₂′ − x₂), x₁ = τ_{X₁}/2h², x₂ = −τ_{X₂}/2h⁴ — a Q-form via x = τ₁(u₀)τ_F-type quantities of `even-bridge-lpairs` §9c (P_K = c₀ + c₁τ_F + c₂xτ_F + c₃y + c₄yτ_F with c₃ = c₄ = 0).
- Open: Q-form of ⟨0_w|Ψ⟩ is still unconnected to any admissible K; a z₂ ≠ z₃ generalisation of K_κ needs ≥ 11 additional localised operators (none known).

## Files

`notebook: reference-admissibility-Kkappa.wb`; `checks/session-2026-10-13/refmap/`: `rm_base.py` (exact pair-level machinery on pgbase), `helpers.py`, `scan.py`/`scan_out.txt`, `obstruction.py`/`obstruction_out.txt`, `genericz.py`/`genericz_out.txt`, `smallint.py`/`smallint_out.txt`, `covariance.py`/`covariance_out.txt`,
`l3check.py`/`l3check_out.txt`, `r6_l3.py`, `kw_on_refs.py`, `symR6.py`, `symR2.py` (sympy symbolic coefficients), `untw.py`/`untw_out.txt`, `twloci_generic.py` (twisted-vacuum loci at generic z).
