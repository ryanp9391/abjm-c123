# even-bridge-lpairs: results (python; Wolfram kernel dead / licence expired)

Scripts here (run `python <script>` in this folder): `eng.py` (fast t1, t2, t3 on bras, diagonal twist, exact mod
P = 268435399 or complex128; validated against pgbase.py exactly), `lp.py` (pair words), `modlin.py`/`modspan.py`
(exact linear algebra mod P), `supp.py` (which pairs a bra differs from W_a on), `task1*.py`.
Point (unless stated): h = 1/2, ABJM point theta_A,n = theta_F,n + 2h with DISTINCT theta_F = (-3/5, -11/8, -7/9, 5/6),
z = (3/2, y, y, 4/7), y = 5/4, a = 1/3.  "exact" = identity mod P (false-positive probability ~1e-8 per entry).

## Task 1 -- l-pair words (all exact mod P)
1. Coordinator's identity confirmed as operator identity (random product bra, L = 2, 3): t2(u0) = t1(u0) t1(thF),
   t3(u0) = 0 (full), u0 = thF + h, ABJM point.  [task1a.py]
2. Pair-1 word reproduced: <W_a^L|K^(1) = (1-k z1)^2 <W_a'(x)W_a^(L-1)| at L = 2, 3, kappa in {3/10, -7/5, 1/z4},
   generic and periodic theta.  [task1a.py]
3. Pair-p word K^(p) (X's at thF_p, lambda over n != p) does NOT localise on pair p for p >= 2: every pair-p word
   (X1, X2, X3, X3X1 at pair p) changes pairs {1..p} (block localisation, as 2603 eq. `localized`).  [task1b.py]
   Naive products prod_{p<=l} K^(p) fail for l = 2 (L = 2, 3) and l = 3 (L = 3).
4. Block localisation (2603-type) holds for all MULTILINEAR monomials prod_p w_p, w_p in {1,X1,X2,X3,X3X1}^(p):
   <W_a^3| m = lambda_m (<W_a^2| m on the 2-pair chain) (x) <W_a|; same-pair squares (X1^2, X1X2, X2^2, X2X3, X3^2)
   do not block-localise.  [task1e.py]
5. TARGET NOT REACHABLE: <W_a'^l (x) W_a^(L-l)| is NOT in the span of the 5^l multilinear monomials (rank 25 / 125,
   full rank, inconsistent) for (L,l) = (2,2), (3,2), (3,3), kappa = 3/10 and 1/z4.  [task1c.py]
   Not in the span of per-pair-degree <= 2 monomials (100, rank 100) either; per-pair degree 3 (400 monomials)
   spans the whole 285-dim cyclic module of <W_a^3| (so "reachable" only trivially, non-local).  [task1i.py]
6. Structure: the W-type product states in the 2-block span (multilinear OR degree <= 2) form the 5-dim space
   span{u_i}(x)W_a  +  W_{a r}(x)span{u_j},  r = z1/z4  (u0 = e1e4, u1 = e4e4 - e1e1, u2 = -e4e1, W_b = u0 + b u1 + b^2 u2):
   pair 1 rotated freely with pair 2 untouched, OR pair 2 rotated freely with pair 1 twisted by the FIXED Lambda^{-1}
   (a -> a z1/z4).  No W_b (x) W_c with both b, c != a, a r.  [task1f.py, task1j.py]
7. Right-end words exist: Y1 = t1(thA_L), Y2 = t2(thA_L), Y3 = t3^(thA_L + h) localise on pair L; {1,Y1,Y2,Y3,Y1Y3}
   localise, rotate pair L freely.  Left x right products: W-type intersection again only
   {W_b (x) ... (x) W_a} + {W_a (x) ... (x) W_c} (dim 5) -- first and last pairs cannot be rotated together. [task1h.py]
8. Pretwisted pair words (exact mod P, L = 3, 4, every p, kappa in {3/10, 1/z4}):
   <W_{ar}^{p-1} (x) W_a^{L-p+1}| K^(p)_kappa = (1-k z1)^2 <W_{ar}^{p-1} (x) W_a' (x) W_a^{L-p}|,  r = z1/z4,
   (X3X1)^(p) = K^(p) at kappa0 = 1/(z1+z4) (s = 0): <..|(X3X1)^(p) = 24 h^7 y^2 l1 l3 z4^2 <..W_{ar} at pair p..|,
   chain word (X3X1)^(1)...(X3X1)^(p-1) K^(p) on <W_a^L| = c (1-kz1)^2 <W_ar^{p-1} W_a' W_a^{L-p}|.
   => prod_p K^(p)_{kappa_p} localises/rotates iff kappa_1 = ... = kappa_{l-1} = 1/(z1+z4) (fixed Lambda^{-1}
   twist on the pairs left of the free one): the alternating-chain loc3 analogue is the END-GROWING chain, not a
   common-kappa block.  L = 4, l = 2: target <W_a'^2 W_a^2| not in the 25 multilinear monomials (rank 25). [task1k.py]

## Task 2 -- even bridges at L = 3 through the untwisting limit [task2.py, untw.py; numeric complex128]
Twist z = 1 + eps phi, phi = (0.9, .35, .35, -1.1), theta_F,n = eps (0, .37, -.52), theta_A = theta_F + 2h; Psi0 = the 26
HW cyclic states of data/abjm-hw-states.wxf (L3) (conventions checked: eigen-residual 1e-14, C reproduced);
twisted state = max-overlap eigenvector in the weight sector (overlap 1.000000 for all 26 at eps = 1e-4, gaps >= 0.038).
Formula <W_a'(x)W_a^2|Psi> = K_kappa(tau) <W_a^3|Psi>/(1-kappa z1)^2: |RHS - LHS| <= 1e-6 (cancellation) everywhere.
Families: F1 kappa = 1/(z1 - a(z1-z4)), a = 1/3 (a'-a = -1) -> C4 = <W_-1 vac^2|Psi0> (M=1);
F2 kappa = 1/(z1 + a(z1-z4)), a = 1/3 (a'-a = +1, shift by -a') -> C2 = <vac W_-1^2|Psi0> (M=2);
F3 kappa = 1/z4, a = -1 -> C2 (M=2);  F4 kappa = 1/z4, a = +1, shift by -1 -> C4 (M=1);
G global rho(1-kappa Lambda) (M = 3): P = [(1-kz4)/(1-kz1)]^3 prod_i (1-kz_i)^{mu_i} -> C0.
ALL converge O(eps) to the data values (errors at eps = 1e-2/1e-3/1e-4, e.g. state 16 [1,0,1] tgt -0.171126:
F1 6.0e-3/6.0e-4/6.1e-5, F2 3.8e-3/3.8e-4/3.8e-5, F3 = F4 1.5e-3/1.5e-4/1.5e-5, G 3.3e-5/3.3e-7/3.3e-9 (tgt 0)).
Nonzero targets at L = 3: vacuum (1) and the six [1,0,1] states (C2 = C4 = -0.171126, -0.027138, -0.035872 -+ 0.028637i,
0.129924 -+ 0.081713i); all other 19 states target 0, reproduced to O(eps).
Symmetry used: for HW Psi0, C_M = (-1)^{L+mu4} C_{L-M} (HW shift by +1 then rho(diag(1,1,1,-1))), so the one-pair word
covers M = 1 and M = L-1 at every L; at L = 3 that is every even bridge.  M = 2 at L = 4 is the first one NOT covered.

## Task 3 -- direction independence and the z2 = z3 constraint
3a [task3a.py, numeric]: kappa = 1 - c eps exactly, a = (phi1 - c)/(phi1 - phi4) (a'-a -> -1 at leading order; derived from
   a'-a = a kappa (z1-z4)/(1-kappa z1) -> a (phi1-phi4)/(c-phi1)); phi in {(.9,.35,.35,-1.1), (-.7,1.3,1.3,.4), (1.5,-.6,-.6,.2)},
   c in {0.3, 2.0}; states 0,4,16,17,18,20,22,25: ALL six (phi, c) runs converge O(eps) to the same C4 values, e.g. state 16
   (tgt -0.171126): |RHS-tgt| at eps=1e-4 = 3.4e-5, 1.9e-4, 3.4e-5, 1.1e-4, 1.8e-5, 3.2e-4; |RHS-LHS| <= 2e-6.
3b [task3b.py]: z2 != z3 at finite generic twist (L = 3, z = (3/2, 5/4 +- d/2, 4/7)): relative residual of the best
   word is LINEAR in d = z2 - z3: K_A (y -> (z2+z3)/2 in the X2 term, y^2 -> z2 z3) 5.25e-3, 5.25e-4, 5.25e-5, 5.25e-6
   for d = 1e-1..1e-4 (kappa = 3/10), = least-squares optimum over all 5 coefficients (5.19e-3, 5.25e-4, ...);
   K_B (y = z2) 6x worse.  => z2 = z3 is a genuine restriction at finite twist (obstruction O(z2 - z3)).

## Task 4 -- Q-form of the pair-p factor (ABJM point; polynomial gauge of sov_l4: Baxter sum_a (-1)^a tau_a(u) Q(u-ah) = 0,
## Q_i = z_i^{u/h} Phi q_i, Phi(u) = V4(u) Phi(u-h), V4(u) = prod_A(u-thA) prod_F(u-thF-h)) [qf.py, qdict.py, task4a-d.py]
u0 = thF + h = thA - h (pair 1).  "rest" = product over the OTHER pairs only; F3(u) = prod_F (u-thF)(u-thF-h);
P_I(u) = det[z_i^{2-s} q_i(u - s h)]_{i in I, s=0,1,2}, q_I = P_I/F3 (reduced Casoratian = dual Q, polynomial), Z_I = prod_{i in I} z_i.
  tau_F = tau1(thF) = Q_j(thF)/Q_j(thF-h) = z_j V4(thF) q_j(thF)/q_j(thF-h)            (any j)
  x = tau1(u0) = [V4rest(u0)/F3rest(thF)] Z_I^4 det[z_i^-1 q_i(thF+h), z_i^-3 q_i(thF-h), z_i^-4 q_i(thF-2h)] / q_I(thF)
  tau2(u0) = x tau_F  (operator identity t2(u0) = t1(u0) t1(thF));  Q-form: same prefactor x V4(thF), columns (1,2,4)
  tau3^(u0) = h V4rest(u0) V4(thF) V4(thF-h) Z_I [F3rest(u0)/F3rest(thF)] q_I(thF+h)/q_I(thF)     (any triple I)
  i.e. X1 = single-Q ratio (F site), X3 = single DUAL-Q ratio q_I(thF+h)/q_I(thF) (A site), X2 = x X1 with x a gapped
  Casoratian / reduced Casoratian.  Derivation: delta -> 0 limit (thA = thF + 2h + delta) of the generic A-site Cramer
  formulas tau_a(thA-h) = D_a/D_0 (the zeros of D_0 and V4(thA-h) are vacuum factors, cancelled explicitly).
  Checks: L = 2, all 256 states x 4 triples, max rel. dev. 5.6e-12 (tau_F), 3.3e-8 (x, tau2, tau3^); L = 3 sectors
  (2,0,0,-2) [24 states] 1e-6, (2,1,-1,-2) [9 states] 2e-8 (float64, Baxter-conditioning limited).  Continuity: at
  delta = 1e-2, 1e-3, 1e-4 tau1(thA-h) -> x and tau3(thA-h)/delta -> tau3^(u0) linearly in delta.
  Factorisation: K(kappa) = 1 - kappa(2S - Sig) + kappa^2 (S^2 - S Sig + Pi) (S = z1+z4) does NOT factor into
  rational linear factors: at L = 2, z2 = z3, 53 of 256 states with real eigenvalues have complex roots.
3c [task3c.py]: untwisting with z2 - z3 = 0.5 eps^q (q = 1, 2, 3; i.e. q = 1 is a direction with phi2 != phi3), word K_A,
   families F1 (tgt C4) and F3 (tgt C2), states 0,4,16,17,18,20,22,25: the finite-twist identity is violated only at
   |RHS-LHS| ~ 0.1 (z2-z3)^2/eps (state #4 [2,0,2]: 2.6e-4, 2.7e-5, 2.7e-6 for q = 1; 2.6e-8, 2.7e-11 for q = 2;
   [1,0,1] states and singlets: at the numerical noise level <= 8e-7) and the LIMIT IS UNCHANGED: |RHS - tgt| identical
   to the z2 = z3 runs (O(eps), e.g. #16: 3.8e-3/3.8e-4/3.8e-5 for every q).  => z2 = z3 is needed for exactness at
   finite twist, but the limit only needs z2 - z3 = O(eps) (any phi; z2 - z3 >> sqrt(eps) would spoil it).

## Task 5 -- summary (also the last cell of Paul/Mathematica/su4-alternating/even-bridge-lpairs.wb, cells unrun)
<W_a'(x)W_a^{L-1}|Psi_eps> = K_kappa(tau) <W_a^L|Psi_eps>/(1-kappa z1)^2 exact (z2 = z3); K in Q's per Task 4.
M = 1: kappa with a'-a = -1 (any a != 0; or kappa = 1 - c eps, a = (phi1-c)/(phi1-phi4)), or kappa = 1/z4, a = +1.
M = L-1: kappa = 1/z4, a = -1, or a'-a = +1 (HW shift by -a').   M = L: global weight factor.  C_M = (-1)^{L+mu4} C_{L-M}.
2 <= M <= L-2 (first case L = 4, M = 2): NOT reachable -- no common-kappa l-pair word (exact rank tests), only end-growing
chains with fixed Lambda^{-1} on the pairs left of the free one.  Remaining: <W_a^L|Psi> in SoV/Q form, middle bridges, norms.
