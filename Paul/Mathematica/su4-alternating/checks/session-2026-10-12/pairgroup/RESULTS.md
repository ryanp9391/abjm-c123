# pair-group-rotation: results (python; Wolfram kernel unavailable)

Scripts in `py/` (run with `python <script>`; `pgbase.py` = fused t_k for k = 0..4, exact Fractions or floats).

## Conventions
Chain A F A F ..., bras = row vectors, L_s(u) = (u - th_s) + h sum e_ij (x) rho_s(E_ij), rho_A(E_ij) = e_ij, rho_F(E_ij) = -e_ji,
T = L_1...L_N G.  t_k(u) = prod_s g_s(u) tr[L^(k)_1 ... L^(k)_N wedge^k G],  L^(k)_s = (u - th~_s) + h sum_ab D_k(e_ab) (x) rho_s(E_ab),
D_k = derivation action on wedge^k C^4, th~_A = th, th~_F = th + (k-1)h, g_A = prod_{j=1}^{k-1}(u-th-jh), g_F = prod_{j=0}^{k-2}(u-th-jh).
Verified against quantum minors for k = 1..4 (random twist, 2 and 4 sites, 1e-15) [`verify_fused.py`].
"hat" t_k = t_k with exactly vanishing scalar factors (u - th - jh) = 0 dropped.

## 1. Derivation (2603 App. K for the alternating chain)  [`pairvac_pattern.py`, `scan_d.py`, exact/numeric]
- <e1 e4|(L_A L_F)_jk(u) = (u-thA)(u-thF) d_jk <e1e4| - h(u-thA) d_k4 <e1 e_j| + h(u-thF) d_j1 <e_k e4| - h^2 d_j1 d_k4 <S|  (exact, symbolic).
  Zero pattern: identically 0 for (2,1),(3,1),(4,1),(3,2),(4,2),(4,3) AND (2,3); extra zeros (2,4),(3,4) at u = thA, (1,2),(1,3) at u = thF.
- Column-support lemma: rows I of the wedge^k auxiliary space with 1 not in I, 4 in I (set D_k) are diagonal on BOTH vacua <e1| (A) and <e4| (F),
  with I-independent eigenvalues; so pairs 2..L in the vacuum are "transparent" on D_k. If pair 1's Lax product at u0 has column support in D_k,
  then <vac^L| t_k(u0) = lambda_k(u0) <vac_pair| t_k^(1 pair)(u0) (x) <vac^(L-1)|  for ANY twist, hence for any global rotation.
  F site at u0 = th~_F kills all columns without 4; A site at u0 = thA - h kills all columns containing 1. Therefore:
  k=1: t1(thF) (any geometry); k=3: t3^(thA - h) (any geometry); k=2: t2(thF + h) ONLY IF thA - thF = 2h (the ABJM point).
  lambda_k(u) = prod_{n>=2} prod_{j=0}^{k-1}(u-thA_n-jh) * (u-thF_n-kh) prod_{j=0}^{k-2}(u-thF_n-jh)   (vanishing factors dropped).
- Scan (4 sites, random twist, random BPS ref), d = (thA-thF)/h in -5..5, t_k(thF+mh), k=1..3, m=-6..7: only t1(thF), t3(thA-h), and at d = 2 only also t2(thF+h).
- Every localising word equals lambda x (same word on the 1-pair chain), exactly [`pairword_check.py`, Fractions, random rational twist & BPS ref].

## 2. Task 1 (good companion reference <0| = e1 e4 e1 e4, 4 sites, exact)  [`task1_scan.py`, `task1_struct.py`]
Single-site localisation, generic theta: t1(th2) -> site 2: hth12th23(th24-h) <0|rho2(G) (= -2263/14000, vector e4 G_c = e3);
t1(th3) -> site 3 by G^T; t3(th1-h) -> site 1 by G^-T; t3(th4+3h) -> site 4 by G^-1. NO t2(th2+s) and NO t3(th2+s) localises on site 2
(t2(th2+h) excites sites {2},{1,2},{2,3}; t3(th2+3h) excites {2},{3},{4},{2,3}). At the ABJM point t2(th2+h) excites only {2},{1,2} (pair 1).
=> The F-site dictionary M2, M3 does not exist at 4 sites; the wedge structure lives on the PAIR, not on the F site.
Same-point generating functions sum_k c_k t_k(u0) (u0 = thF+mh, thA+mh, m=-3..4, also "staircase" u0+(k-1)h): no combination localises
beyond the individually localising terms (generic theta and ABJM) [`samepoint.py`].

## 3. Localised algebra on pair 1 at the ABJM point  [`span_full.py`, `pair_sym.py`, `pair_rank1*.py`]
All individually localising words (products up to length 4 of t_k(thF+mh), k=1..3, m=-2..5): span = 5 dims on any reference:
{1, X1 = t1(thF), X2 = t2(thF+h), X3 = t3^(thF+h) (= t3^(thA-h)), X3 X1}.  (generic theta: 4 dims, no X2.)
On the pair vacuum with general twist G (matrix M_ij, A index i, F index j):
X1 = 2h^2 e1 (x) e4 G - h^2 G41 S ;  X3 = 12h^5 (e1 cof G) (x) e4 ;  X3X1 ~ rho_pair(G^-T) ;  X2 = 2x2 minors of G with rows {2,4},{3,4} (no group form).
Rank-one (product) elements of the span: only v(a + b G^-T) (x) u and v G^-T (x) u(a + b G): NO common group rotation, for any reference with u.G.v != 0.
Eigenvalues on weight states e_i (x) e_j (diagonal twist): x1 = X1/(2h^2) -> z_j, x2 = -X2/(2h^4) -> z_j (e1(z) - z_i - z_j), x3 = X3/(12h^5) -> det/z_i.
rho_pair(1 - kappa Lambda) needs (1-kappa z_i)/(1-kappa z_j) on weight states and 1 on zero-weight states.

## 4. THE TEST on non-good references (diagonal frame)
- Generic BPS (v (x) u)^2 [z2 != z3 or z2 = z3], ABJM point: localised span 5, span+target 6 (exact) -> no common rotation.
- W_a^2, ABJM point, z2 != z3: span 5, span+target 6 (exact). Obstruction = zero-weight components E22, E33 of X1 (prop. E22+E33) and
  X2 (prop. z3 E22 + z2 E33): the cyclic space of W_a on the pair is 6-dim (2 weight states + 4 zero-weight eigenstates), the algebra gives 5.
- W_a^2, generic theta, z2 = z3: span 4, span+target 5 (exact).
- **W_a^L, ABJM point, z2 = z3 = y: WORKS.**  The zero-weight eigenstate odd under 2<->3 decouples, cyclic space 5 = algebra 5.

### The word (exact, closed form, a-independent; pair level solved symbolically in z1, y, z4, h, a, kappa) [`pair_solve.py`]
  K_kappa = s^2 + kappa s ( X1/(h^2 l1) + X2/(2h^4 y l2) + X3/(6h^5 y^2 l3) ) + kappa^2 X3 X1/(24 h^7 y^2 l1 l3),   s = 1 - kappa(z1+z4),
  X1 = t1(th2), X2 = t2(th2+h), X3 = t3^(th2+h), l_k = lambda_k (formula above; L-dependent only through pairs 2..L).
  <W_a^L| K_kappa = (1 - kappa z1)^2 <W_a' (x) W_a^(L-1)|,   a' = a (1-kappa z4)/(1-kappa z1)   (pair-level: K = (1-k z1)(1-k z4) rho_pair(1-k Lambda) on W-pairs).
  Normalised: K = s^2 + kappa s (2x1 - x2/y + 2x3/y^2) + kappa^2 x3 x1/y^2.
  Exact verification (Fractions): L = 2 and L = 3 pairs, a in {1/3, -1}, kappa in {3/10, -7/5, 1/z4} (incl. a' = 0) [`word_verify.py`]. Locality holds.
  kappa = 1/z4 (vacuum pair):  <vac (x) W_a^(L-1)|Psi> = (z1^2 - z1 Sigma + Pi)/(z1 - z4)^2 <W_a^L|Psi>, Sigma = 2x1 - x2/y + 2x3/y^2, Pi = x3 x1/y^2 (eigenvalues).
- Same K on a generic BPS reference (z2 = z3): pair state rank 2 (entangled), traceless, 3.5% from d1 d4 rho(D) [`negcontrols.py`].

### Combinations of NON-localising words [`combo_search.py`, `combo_L3.py`, `sparse_generic.py`]
At L = 2 universal (3 values of a) combinations of 253 words reach W_a' (x) W_a even for z2 != z3 and generic theta (the 2-pair algebra is too rich).
At L = 3: ABJM & z2 != z3: rank 5, target residual 2.3e-2 (fails); generic theta & z2 != z3: rank 4, fails; ABJM & z2 = z3: reaches (5);
generic theta & z2 = z3: reaches, but the L1-sparsest combination has 80 words with coefficients 1e-3..1e7 (no structure; L-universality untested).

## 5. Untwisting (numeric, complex128) [`untwist.py 2`, `untwist.py 3`]
Homogeneous ABJM chain thF,n = eps*thdir_n, thA = thF + 2h (thdir = (0, .37, -.52)), z = 1 + eps*(0.9, 0.35, 0.35, -1.1) (z2 = z3),
untwisted HW cyclic joint t1,t2,t3 eigenstates (L=2: 6, L=3: 26) computed in python; twisted state = eigenvector in the same weight sector
with max overlap (all non-degenerate, overlap -> 1).
(i) kappa = 1/z4, a = -1: RHS = K(tau)<W_-1^L|Psi>/(1-z1/z4)^2 = <vac (x) W_-1^(L-1)|Psi> (|RHS-LHS| <= 1e-6 machine cancellation; K/N ~ 1e8, <W|Psi> ~ 1e-10).
(ii) kappa = 1/(z1 - a(z1-z4)) -> 1 (a = 1/3, a'-a = -1 exactly), RHS -> <W_-2/3 (x) W_1/3^(L-1)|Psi0> = <vac^(L-1) (x) W_-1|Psi0> (HW trick + cyclicity).
Errors |RHS - target| at eps = 1e-2, 1e-3, 1e-4 (all states, O(eps)):
 L=2 (i): singlets -1/sqrt40, 1/sqrt24: 1.6e-3/1.6e-4/1.6e-5, 2.0e-3/2.0e-4/2.0e-5; [1,0,1]: 5.4e-3/5.4e-4/5.4e-5 (target 0); vacuum 1 (2.6e-13).
 L=2 (ii): 5.2e-4/5.3e-5/5.3e-6 ; 7.0e-4/6.8e-5/6.8e-6 ; [1,0,1] 3.1e-3/3.1e-4/3.1e-5.
 L=3 (i) and (ii): targets (C2 = C4 here) -0.153484 (x2), -0.035872 -+ 0.028637 i, -0.171126, -0.027138 for the [1,0,1] states, 0 otherwise;
 max error 3.8e-3 / 3.8e-4 / 3.8e-5.
