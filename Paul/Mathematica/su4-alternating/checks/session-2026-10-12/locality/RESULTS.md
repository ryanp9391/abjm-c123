# Pair-rotation locality test at L = 3 pairs: results (python; Wolfram kernel unavailable)

Setup exactly as specified (chain A F A F A F; bras as row vectors; Lax L_A = (u-th)+hK, L_F = (u-th)-hP; diagonal twist Lambda = diag(z);
t1, t2 (6-dim wedge^2 aux, F shifted by h), t3 (dual Lax, A -h, F +3h, twist det G G^-T) as in dual-lax-t3.wb; h = 1/2,
theta = (2/7,-3/5,8/11,-11/8,5/13,-7/9), z = (3/2,5/4,2/3,4/7), reference psi0 = W_{1/3}^{(x)L}).
Sanity (done): [t_a,t_b] = 0 on random vectors (4 and 6 sites, ~4e-16 relative, floats; exact mod p in the modular code);
calibration psi0 = c1 (psi0 rho_1(Lambda^-1)) t1(theta_1) holds EXACTLY (Fractions) for psi0 = W_a^{(x)L} at 4 and 6 sites, three theta points
(NOTE: it is false for a general product bra or for W_{a'} (x) W_a^{L-1}; "any product bra" in the task text must be read as "any globally
rotated vacuum e^{a E14} vac^L").  Speed: 4 ms per t3 application on a 4096-dim vector in numpy (task benchmark 0.2 s was for WL).

## 1. Word lists

* factors t_a(theta_s + m h), a=1,2,3, s=1 (A site of pair 1), 2 (F site of pair 1), m=-4..4: 54; 6 vanish identically (zeros of the prefactor g:
  t2(th_A+h), t2(th_F), t3(th_A+h), t3(th_A+2h), t3(th_F), t3(th_F+h)) -> **48 non-zero factors** (24 per site). Same list at 4 and 6 sites.
* words <= 2 factors (1 + 48 + 1176 incl. squares) = **1225**; triples A F F: 24 x 300 = 7200; A A F: 7200 (task classes: 15625 words total);
  A A A and F F F triples: 2600 each (all-type <=3 factors: 20825 words); (1, t1(th1)^{+-1}) x t1(th2+m1 h) t2(th2+m2 h) t3(th2+m3 h): 504 x 3 = 1512 words.
* inverses: t1(th1)^-1 psi0 = c1 psi0 rho_1(Lambda^-1) is explicit (exact); for every other inverse, inverse words were handled inside the
  exact cyclic space (inverse of a restricted operator is a polynomial in it).

## 2. Spans, localised subspaces, targets (EXACT: arithmetic mod p = 4194301, confirmed with p = 4189277; ranks are exact up to a ~1e-4 accident probability, none observed between the two primes)

Layered span V_l = span of all words with <= l factors (all types, positive powers, 48 factors). "loc" = dim of V_l intersect (C^16 (x) W_a^{L-1}) = bras of the form x (x) w0 on pair 1.

| L (sites) | l=0 | l=1 | l=2 | l=3 | l=4 | l=5 | l=6 | full algebra orbit (Krylov) |
|---|---|---|---|---|---|---|---|---|
| 2 pairs (4 sites) rank | 1 | 11 | 41 | 53 | 54 | 54 | | **54** |
| loc | 1 | 3 | 6 | 6 | 6 | 6 | | **6** |
| T1,T2 (all 4+3 values) in V_l | no | no | **yes** | yes | yes | yes | | yes |
| 3 pairs (6 sites) rank | 1 | 17 | 117 | 400 | 563 | 587 | 588 | **588** |
| loc | 1 | 3 | 4 | 4 | 6 | 6 | 6 | **6** |
| T1,T2 in V_l | no | no | no | **no** | **yes** | yes | yes | yes |
| T3 controls in V_l (T3a / T3b) | no / no | yes / no | yes / yes | yes / yes | | | | |

Task word classes at L = 3 (2 pairs words + AFF/AAF triples = 15625 words): rank 400, loc 4, T1/T2 absent (exactly); adding the 504 F-site t1 t2 t3 products times {1, t1(th1), t1(th1)^-1}: see section 7 (rank unchanged/changed as listed there).
The result is the same at three theta points (standard; theta' = (1/3,-2/5,3/4,-5/4,2/9,-7/10), which is NON-generic because theta_3 - theta_4 = 4h
(cyclic dims 48 / 534 instead of the generic 54 / 588, loc 6, layer structure l = 2 / l = 4 unchanged: ranks 1,10,36,47,48 and 1,16,106,359,510,533,534);
and a generic third point (19/23,34/29,41/40,-19/12,-7/11,13/17): L=2 rank/loc as standard, L=3 ranks 1,17,117,400,563,587,588 identical).

Float (double) cross-check of the same quantities: ranks 117 (l<=2), 400 (l<=3) agree with the exact ones when the singular-value gap is read correctly, but
the float word matrix is badly scaled (column norms 1e-8 .. 1e16, noise columns) and naive thresholds (1e-9 on normalised columns) give wrong ranks (411, 440).
So all rank / membership claims above are the modular ones; float is used only for exploration. mpmath/Fraction spot checks: section 6.

## 3. The structure behind it (exact identities, Fractions, 4 and 6 sites, 3 theta points)

* psi0 t1(th1)^-1 = c1 psi0 rho_1(Lambda^-1), 1/c1 = h prod_{F s} th_{1s} prod_{A s>=2}(th_{1s}+h).
* psi0 t3(th1-h) = c3 psi0 rho_1(Lambda^-1), c3 = -(3/4) det(Lambda) prod_{F s} (th_{1s}-h)(th_{1s}-2h)(th_{1s}-4h) prod_{A s>=3}(th_{1s}-h)(th_{1s}-2h)(th_{1s}-3h)   (h = 1/2; -3/4 = -6h^3).
* psi0 t1(th1)^-1 t1(th2) = c12 psi0 rho_1(Lambda^-1) rho_2(Lambda), c12 = prod_{s>=3} [A: th_{2s}/(th_{1s}+h); F: (th_{2s}-h)/th_{1s}].
  All three constants are products of per-site (per-pair) factors: each added pair multiplies them by explicit theta-factors: universal / local.
* psi0 t1(th2) = alpha [ psi0 rho_2(Lambda) + (sigma/alpha) S (x) w0 ],  S = sum_i e_i (x) e_i (pair-1 singlet), w0 = W_a^{L-1},
  alpha = h th12 th23 (th24-h) [(L=2; times prod_{s>=5} per pair)], **sigma/alpha = h a (z1 - z4)/th12, independent of L** (checked L=2,3).
  So the F-site step always drags a singlet.
* The 8 individually localised words (<= 3 factors) are IDENTICAL at L=2 and L=3: 1, t3(th1-h), t1(th2), t1(th1)t3(th1-h), t3(th1-h)t1(th2),
  t1(th2)t3(th2+3h), t3(th1-h)t1(th2)t3(th2+3h), t1(th1)t3(th1-h)t1(th2); their x-vectors span only 4 dimensions:
  psi0, rho_1(Lambda^-1)psi0, rho_1(Lambda^-1)rho_2(Lambda)psi0, and rho_2(Lambda)psi0 + (sigma/alpha) S (x) w0.
* The 10 two-word linear relations between words (parallel bras, e.g. t1(th1-h) ~ t2(th1)t3(th1-h), t2(th1) ~ t1(th1-h)t1(th1), t3(th1) ~ t1(th1)t2(th1-h),
  t2(th1-h) ~ t3(th1-h)t3(th1), t1(th1-h)t3(th1) ~ t2(th1-h)t2(th1), and the same five with th2 -> th2+h or +2h pattern) are fusion relations; they hold at L=2 and L=3 identically.
* EXACT (mod p, 3 theta points, L=2,3): T1(a') and T2(kappa) lie in span{psi0, psi0 t3(th1-h), psi0 t3(th1-h)t1(th2), psi0 t1(th2)} (+) C S(x)w0, and S(x)w0 is not in the span of the 4 words.
  **The whole obstruction is the isolated pair singlet S (x) w0** (equivalently the pure F-site rotation psi0 rho_2(Lambda)).

## 4. When is S (x) w0 / the pure F-site rotation reachable?

| object | L=2 | L=3 |
|---|---|---|
| psi0 rho_1(Lambda^{+-1}), v(1-kappa Lambda^-1) x u (T3a) | l = 1 | l = 1 |
| rho_1(Lambda^-1) rho_2(Lambda) psi0, vLambda^-1 x u(alpha+beta Lambda) (T3b) | l = 2 (single word t3(th1-h)t1(th2)) | l = 2 (same word) |
| pure F-site rotations v x u Lambda^k (k = -1,1,2,3), singlet S x w0, rho_1(Lambda^{+-1}) S x w0, T1, T2 | **l = 2** | **l = 4** (not l = 3) |

* Needed degree: l = 2 at L = 2 pairs, l = 4 at L = 3 pairs (two data points; pattern 2(L-1), not tested at L = 4).
* A-site factors alone (24) or F-site factors alone (24) already suffice (same layer thresholds); with m in {-1,0} only (10 factors) the threshold is still 2 / 4;
  t1 only or t3 only never reach T1/T2; t2 only needs l = 4 / 6; (t1,t3) l = 4 / 6, (t1,t2) and (t2,t3) l = 3 / 5.

## 5. Sparse / minimal word sets (L = 2 pairs, the only length where T1, T2 are reachable by few-factor words)

* Residual in the span of all 1225 words (rank 41): 1e-15 for every T1, T2 target; singles only (rank 11): 1.7e-2 ... 6.3e-1.
* L1 (HiGHS, unit-norm columns): support = 41 = rank (no sparsity); iteratively re-weighted L1: still 41; orthogonal matching pursuit: no exact solution with <= 12 words
  (residual 1.6e-1 .. 4e-2 after 8 words); backward elimination from the full pool, 60 random orders: irredundant supports of size 41-42 for T1 and T2 (always a full basis of V_2).
  Exhaustive subset search over {1, 48 singles}: no solution (singles span only 11 dims and do not contain T1/T2).
* The pieces: rho_1(Lambda^-1) rho_2(Lambda) psi0 is a ONE-word target (t3(th1-h) t1(th2)); the pure F-rotation and the singlet have L1 support 41 and OMP residual
  0.08 / 0.15 after 8 words. So T1/T2 = 3 universal local words + (dense, many-word) singlet isolation.
* (MIP: see section 8.)

## 6. Precision control

* Exact rational arithmetic (python Fractions) at 4 sites (all 48 singles and the 41 independent <=2-factor words built exactly) + 30-digit modified Gram-Schmidt (mpmath, dps = 30):
  relative residual of T1/T2 in span(<=1 factor words, 11 independent) = 0.63, 0.12, 0.090, 0.31 (T1) and 0.049, 0.017, 0.13 (T2): clear failure, not a rounding artefact;
  in span(<=2 factor words, 41 independent) = 1e-23 ... 7e-23 (success, limited only by the conditioning of the 41 words); T3a: 5e-30 at l = 1.
* Exact Fractions at 6 sites (single passes, 15-25 s each): calibration, c12, c3 identities of section 3 verified exactly at theta standard, theta', generic theta''.
* Floats are not used for any claim of non-membership at 6 sites (conditioning); those are the exact mod-p rank statements (two primes).

## 7. Inverse-word and 4-factor classes (exact, mod p, in the cyclic space; L = 3 unless stated)

| word set | L=2: rank / loc / T1,T2 | L=3: rank / loc / T1,T2 |
|---|---|---|
| 1 + 48 positive singles | 11 / 3 / no | 17 / 3 / no |
| 1 + 48 INVERSE singles t_f^-1 (all 48, exact inverse inside the cyclic space) | 49 / 3 / no | 49 / 3 / no |
| 1 + positive singles + inverse singles (97 words) | **54 / 6 / yes** | **61 / 3 / no** |
| ... + (inverse f)(single g) (2401 words) | 54 / 6 / yes | **588 / 6 / yes** |
| 1 + singles + {t1(th1)^-1, t1(th2)^-1} x (1,single) [task inverses] | 24 / 4 / no | 42 / 4 / no |
| <=2 positive + {t1(th1)^-1, t1(th2)^-1} x (<=2) | 51 / 6 / yes (L=2) | 229 / 4 / no |
| <=3 positive + {t1(th1)^-1, t1(th2)^-1} x (<=2) | 53 / 6 / yes | 400 / 4 / no |
| task classes (1225 + AFF + AAF = 15625 words) | 53 / 6 / yes | 400 / 4 / no |
| ... + 504 FFF t1 t2 t3 products x {1} | 53 / 6 / yes | 400 / 4 / no |
| ... + x t1(th1), x t1(th1)^-1, both (4-factor class) | 54 / 6 / yes | 454, 454, 506 / 4 / **no** |
| <=2 positive + all 48 inverses x (1,single) (3577 words) | 54 / 6 / yes | **588 / 6 / yes** |
| inverses of t1 only x (1,single), t2 only, t3 only, A-site only, F-site only, m=0 only (L=3) | | 296 / 4 / no, 269 / 3 / no, 232 / 4 / no, 385 / 4 / no, 385 / 4 / no, 72 / 4 / no |

So **with inverses too the word length grows with L**: at L = 2 single-factor words (positive and inverse, 97 words, rank 54) already reach T1, T2; at L = 3 they do not (rank 61, loc 3); words of length 2 with one inverse factor (2401 words) reach everything (588, loc 6).
Caveat for the inverse words: t_f^-1 psi0 for generic f is NOT local (e.g. t1(th2)^-1 psi0 = -h^-1 psi0 R21^-1 G2^-1 R24^-1 R23^-1: R24^-1 = (th - hP)/(th^2-h^2) acts on psi0 rho2(Lambda^-1) and moves the F-site vector to site 4); only t1(th1)^-1 (c1) is explicit. A float run of the L = 2 all-inverse span gave rank 50-52 instead of the exact 54 and residuals 1e-3..1e-2 for T1, T2: floats are unreliable for inverse words.

## 8. Exhaustive small-support search (float screening, pinv, threshold 1e-7), L = 2

* pool = 105 words (|m| <= 1, <= 2 factors; rank 41): NO subset of <= 4 words contains T1(2/3), T2(1/3), the pure F-site rotation v(x)u Lambda(x)w0 or the singlet S(x)w0 (105, 5460, 187460, 4.6e6 subsets). k = 5 (96e6 subsets, 842 s) also finishes with NO solution for T1(2/3), T2(1/3), the pure F-site rotation or the singlet: **no <= 5-word solution in the 105-word family** (the notebook was written before this finished and says so).
* pool = 97 words {1, 48 singles, 48 inverse singles} (float rank 52 < exact 54): no subset of <= 4 words contains any T1/T2 target (screen only; float conditioning limits this).
* The mixed-integer attempt crashed (memory), so no MIP statement.

## 9. Verdict

At 3 pairs the even-bridge target (T1 / T2) is **not reachable by few pair-1 transfer-matrix words**:
* not by any word of <= 3 factors from the 48 pair-1 factors (rank 400, loc 4; exact, two primes, three theta points), nor by adding t1(th1)^-1, t1(th2)^-1 (<= 2 other factors), nor the (1, t1(th1)^{+-1}) x F-site t1 t2 t3 four-factor class (rank 506, loc 4, T1/T2 absent);
* not by single words incl. all 48 inverses (97 words, rank 61, loc 3), although at L = 2 these 97 words (and the <= 2-factor positive words) do reach T1, T2;
* first reachable at word length 4 (positive; 563-dim span) or length 2 with one inverse (2401 words; 588): the needed length grows from L = 2 to L = 3 (positive: 2 -> 4; with inverses: 1 -> 2); at L = 2 already no <= 4-word solution exists in the 105-word family, L1/OMP/backward elimination give full-rank supports.
* What IS local/universal (same words, same relations, coefficients = products of explicit per-pair theta-factors): the four SoV-type words {1, t3(th1-h), t3(th1-h)t1(th2), t1(th2)} and the T3 controls. T1/T2 reduce EXACTLY to these four words plus the isolated pair singlet S(x)w0 (coefficients of psi0 and S(x)w0 in T1(a'): 9(8a'-7)^2/169 and 35(3a'-1)(8a'-7)/806, identical at L = 2, 3); the singlet is not locally isolable. psi0 t1(th2) = alpha[psi0 rho2(Lambda) + (sigma/alpha) S(x)w0] with sigma/alpha = h a (z1-z4)/th12, independent of L.
* The needed length 2(L-1) (positive words) rests on L = 2, 3 only (L = 4 not attempted).

## 10. Files

Scripts in this folder (python; numpy/scipy/mpmath): lib.py (t1,t2,t3 vector passes, Fraction/float), modlib.py (mod p passes + incremental modular basis), modrun.py / layers_mod.py / krylov_mod.py / cyc.py / cyc_save.py / an_cyc.py / fam_cyc.py / extra_cyc.py / inv_cyc.py / four_plus_S.py (exact modular spans, layers, cyclic space, families, controls, inverse words), exact_ctrl.py / c3sym.py / sigsym.py / coef_exact.py (Fraction/sympy identities and coefficients), mpcheck4.py (mpmath 30-digit), loc_words*.py / circ*.py (localised words, fusion relations), sp*.py / backelim.py / exh*.py / mip*.py (sparse searches), nbcells/ (Wolfram port, UNRUN) and nb_cells.json.
