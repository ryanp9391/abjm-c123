# even-bridge SoV formula -- RESULTS (work in progress, python only; Wolfram kernel dead / licence expired)

Scripts in this folder (copied to the repo: `su4-alternating/checks/session-2026-10-12/sovformula/`).
Base modules copied from this session: `pgbase.py` (pairgroup, fused t_k, bra tensors), `sov_l4.py` (SoV dictionary L = 4 sites),
`core.py`, `exact.py` (twoblock).

## 0. Analytic facts established so far (derivations, checked by hand)
* Diagonal-frame image of the companion SoV vacuum (R = rho(V^{-T}) globally; A: e1 V^{-T} = (lambda_1..lambda_4),
  F: e4 V = (1,1,1,1)):  <0|_diag = (lambda (x) 1)^{(x)L} = lambda_1^L <vac^L| rho(n'D'),  lambda_j = 1/prod_{m!=j}(z_j - z_m),
  r = lambda/lambda_1, n' = ((1,1,1,r4),(0,1,0,-r2),(0,0,1,-r3),(0,0,0,1)), D' = diag(1,r2,r3,1)  (uses sum_j r_j = 0).
  Holds at every L (product state, site by site).
* Ket SoV vacuum |0bar> = e4 e1 e4 e1 (companion) -> diag frame |0bar>_diag = ((1,1,1,1) (x) lambda)^{(x)L}.
* Weight identity: for any weight-w ket Psi (A e_i: +eps_i, F e_j: -eps_j),  <0|Psi> = lambda^w <0bar^T|Psi>,  lambda^w = prod_j lambda_j^{w_j};
  for any weight-w bra PsiL: <PsiL|0^T> = lambda^w <PsiL|0bar>.  (prod_A lambda_{i_A} = lambda^w prod_F lambda_{j_F} on every basis state.)
  => in SoV normalisation (<0|Psi> = 1 = <PsiL|0bar>):  <0|Psi><PsiL|0^T>/<PsiL|Psi> = lambda^w * N_SoV,  1/N_SoV = sum tau~_y M_yx tau_x.
  => if PsiL -> c Psi^dagger (untwisting limit, HW weight space 1-dim): |C_phys(B)|^2 = |P_kappa|^2 |lambda^w| |N_SoV|  (to be verified).
* Coordinator input (ABJM-point agent): at u0 = thA - h = thF + h: t3(u0) = 0 and t2(u0) = t1(u0) t1(thF) exactly, so
  X2 = t1(u0) X1, X3 ~ Z = t3'(u0); regularised basis per pair {1, t1(u0), W, Z} x t1(thF)^k.

## 1. Task 1: localised words on <0|_comp (EXACT, Fractions, t01_pairspan.py; L = 2 pairs; ABJM point pair by pair)
At thA_n = thF_n + 2h, u0 = thF1 + h, <0|_comp = (e1 e4)^L, companion twist G_c(z) (z generic AND z2 = z3 give identical output):
   <0| X = lambda_X * (pair-1 state) (x) (e1 e4)^{L-1}, with pair-1 states (A index, F index), Z-INDEPENDENT:
     1     : e1 x e4
     X1    = t1(thF1)        : 2h^2  e1 x e3
     X2    = t2(thF1+h)      : 2h^4 (e1 x e2 + e2 x e3)        [= t1(u0) X1 exactly on <0|, coordinator's identity]
     X3    = t3^(thF1+h)     : -12h^5 e4 x e4
     X3X1                    : -24h^7 e4 x e3
   lambda_X as in pgbase (lambda_1(thF1), lambda_2(u0), lambda_3(u0), lambda_3 lambda_1).  t1(u0) alone does NOT localise.
* Rank-one (product) elements of this 5-dim span: exactly (a e1 + b e4) x (c e3 + d e4) = <0|(a' + b' X3)(c' + d' X1)
  (the X2 direction is never rank one: its A row e2 forces the coefficient to vanish).  BPS (v.u = 0) iff b d = 0:
  F-only family e1 x (c e3 + d e4) = <0| rho(1 - (c/d) E34), or A-family (a e1 + b e4) x e3.
* K_kappa (pairgroup word) on <0|: pair state  s^2 e1e4 + kappa s [2 e1e3 + (e1e2 + e2e3)/y - 2 e4e4/y^2] - kappa^2 e4e3/y^2,
  Schmidt rank 3 (3x3 minor = -2 kappa^3 s^3 / y^4 != 0)  =>  K_kappa is NOT a common rotation (not even a product) on <0|.
* Physical frame (A: v_phys = v_c V^{-T} (n'D')^{-1} / lambda_1-normalised; F: u_phys = u_c V (n'D')^T): companion F vector
  with polynomial f (e1 <-> x^3, ..., e4 <-> 1) maps to (f_3/lambda_1, r2(f(z2)-f(z4)), r3(f(z3)-f(z4)), f(z4)); companion A
  vector <-> functional l maps to (lambda_1 l(pi_1), lambda_1(l(pi_2)-l(pi_1)), lambda_1(l(pi_3)-l(pi_1)), l(1)), pi_j = prod_{m!=j}(x-z_m).
  Since no localised word puts e1 on the F site, (u3)_1 = p23 = 0 IDENTICALLY for every bra in the localised span:
  x = p23/(a2 d3) = 0, so the singlet (even-bridge, F/F_vac = c x y) structure constants are invisible from <0| with
  pair-1-localised words.  The A-family gives y = p32/(a3 d2) finite only for a/b = c4/z1 + O(lambda_1^{-1}) (odd bridge).

* t02_pairwords.py (exact, L = 2, z2 = z3 and generic z): among the 16 regularised pair-1 SoV words
  <0| O1 t1(thF1)^k (O1 in {1, t1(u0), W, Z}) the localising combinations form EXACTLY the 5-dim span above, and each localising
  combination is a SINGLE regularised SoV word: <0|1, <0|t1(thF), <0|t1(u0)t1(thF) (= X2), <0|Z (= X3), <0|Z t1(thF) (= X3X1).
  W and t1(thF)^{2,3} never localise.  (L = 3 run: t02_L3.log.)  So the SoV expansion of every localised pair-1 bra on <0|
  has at most 5 words, and P_K(tau) = c0 + c1 tauF + c2 x tauF + c3 y + c4 y tauF  (x = tau1(u0), y = tau3'(u0), tauF = tau1(thF1)).

## 2. Norm identity and SoV-normalised formula (VERIFIED, mpmath 80 digits, L = 2, t08_mp_L2.py, t08b_mp_L2.py)
Untwisting: thF_n = eps thdir_n, thA_n = thF_n + 2h (h = 1/2), z = 1 + eps phi, diag frame; saved HW states (data/abjm-hw-states.wxf).
  N_SoV := <0|Psi><PsiL|0bar>/<PsiL|Psi>  (= 1/sum_{x,y} tau~_y M_yx tau_x whenever bra and ket SoV bases are complete).
  herm := |<0|Psi>|^2 / (|lambda^w N_SoV| ||Psi||^2) = 1 + O(eps) for ALL six L = 2 states, three phi, eps = 1e-3/1e-5/1e-7:
     e.g. eps = 1e-7: vacuum 1, [0,2,0] 1.000000296, [1,0,1] 1.000000222 (x2), singlets 1.00000037, 1.000000074.
  => PsiL -> c Psi^dagger in the untwisting limit, and for ANY word K (any polynomial in transfer matrices):
       |<0|K|Psi>| / ||Psi|| = |P_K(tau)| sqrt|lambda^w N_SoV|,   lambda^w = prod_j lambda_j^{w_j}
       |C_K(Psi)| := |<0|K|Psi>/<0|K|vac>| = |P_K(tau)/P_K(tau_vac)| sqrt|lambda^{w - w_vac} N_SoV(Psi)/N_SoV(vac)|,
     N_SoV(vac) = (lambda_1 lambda_4)^L, lambda^{w_vac} = (lambda_1/lambda_4)^L.  Checked: |C_SoVformula| = |C_direct| (1 + O(eps)).

## 3. Physical identification of <0|-based limits: NOT the naive HW reduction
* Naive argument: <B|Psi_eps> = lambda_1^L D'^w <Btilde| rho(n') |Psi_eps>, rho(n') Psi_0 = Psi_0 => C = <Btilde|Psi_0>/<Btilde|vac>.
  FAILS because the physical pair-1 vectors of localised bras have components of order eps^{-2} (A site e2, e3:
  lambda_1 c4 (1/z_{2,3} - 1/z1)), which multiply the O(eps) non-HW corrections rho(n')Psi_1 etc.
* A-site (odd-bridge) word K = alpha + beta X3, alpha = c4/z1 + 1/(lambda_1 y0), beta = -1/(12 h^5 lambda_3)
  (companion pair-1 A vector alpha e1 + e4; naive physical polarisation v3 = (1/y0, O(eps^-2), O(eps^-2), 1), u3 = e4, x = 0, y = y0).
  Naive prediction for [1,0,1]: C = y0 <e4e4e1e4|Psi0> = y0 (0.144338 -+ 0.25 i).  ACTUAL limit (mp, eps -> 0, rate O(eps)):
     C = (y0/3) <e2e2e1e4|Psi0> = (y0/3)(-0.144338 -+ 0.25 i) = (y0/3) e^{-+ i pi/3} c_y,   |C| = |c_y y0|/3,
  independent of phi (three directions) and of thdir ([0,0.37], [1/3,-4/5]), linear in y0 (0.7, -1.3).  Singlets -> 0
  (x = 0), [0,2,0] diverges (mu_v ~ eps^-2, as expected).  => the <0|-based formula computes a well-defined, direction-
  independent limit, but its relation to the physical structure constant is a nontrivial (here 1/3 e^{-+i pi/3}) factor
  coming from higher-order twist perturbation theory; only checked for one state type.

## 4. Q-form of the localised-word eigenvalues (VERIFIED, t16_qform.py, mpmath 50 digits, 41 joint eigenstates,
##    L = 2, ABJM point pair by pair thF = (-3/5, -11/8), h = 1/2, generic diagonal twist z = (3/2, 5/4, 2/3, 4/7))
Baxter (sov_l4 convention) sum_a (-1)^a tau_a(u) Q(u - a h) = 0, Q_i = z_i^{u/h} Phi(u) q_i(u), Phi(u) = V4(u) Phi(u-h),
V4(u) = prod_A (u-thA) prod_F (u-thF-h).  At u0 = thF1 + h = thA1 - h:  Q_i(u0) = 0 (Phi), tau3(u0) = tau4(u0) = tau4'(u0) = 0,
Q_i(u0-h) = tauF Q_i(u0-2h), tau2(u0) = x tauF.  The u-derivative of the Baxter equation gives, for every i = 1..4 (rank 3):
     Q_i'(u0) - x D_i + w Q_i(u0-2h) - y Q_i(u0-3h) = 0,     D_i := Q_i'(u0-h) - tauF Q_i'(u0-2h)
  =>  x = det[Q'(u0), Q(u0-2h), Q(u0-3h)] / det[D, Q(u0-2h), Q(u0-3h)]
      y = det[D, Q(u0-2h), Q'(u0)]        / det[D, Q(u0-2h), Q(u0-3h)]
      w = -det[D, Q'(u0), Q(u0-3h)]       / det[D, Q(u0-2h), Q(u0-3h)]     (rows: any three of Q_1..Q_4)
  tauF = tau1(thF1) = Q_j(thF1)/Q_j(thF1 - h) = z_j V4(thF1) q_j(thF1)/q_j(thF1-h)   (any j)
  x = tau1(u0) (eigenvalue of t1(u0)), y = tau3'(u0) (eigenvalue of Z = t3hat(u0)), w = tau2'(u0) - tau1'(u0) tauF.
  Verified (in the equivalent q-gauge form A_i w + B_i x + C_i y + K_i = 0): max rel. error tauF 4.5e-49, x 2.9e-27,
  w 2.8e-27, y 6.2e-27 (numerical u-derivatives), 4th (redundant) equation residual <= 2e-30.
  => for the localised words on <0|:  P_K = c0 + c1 tauF + c2 x tauF + c3 y + c4 y tauF  is a ratio of 3x3 Q-Casoratians
     with derivative columns (plus single-Q ratios).

## 5. Even with a common rotation the even bridge would be missed; 2603-type dressed references reach it (t17_frames.py, exact sympy)
* The task's expected companion common rotation D_c = 1 - kappa Gc^T is diag 1 - kappa Lam, physical n'(1 - kappa Lam)n'^{-1}.
  Its (1,4) entry is sum_j r_j (1 - kappa z_j) = 0 EXACTLY (sum_j lambda_j z_j^k = 0 for k <= 2; more generally
  (n' f(Lam) n'^{-1})_{14} = f[z1,z2,z3,z4]/lambda_1, the third divided difference).  So the physical A vector has (v3)_4 = 0:
  p32 = y = 0 identically, while the F side gives x = (u3)_1/(u3)_4 = kappa^3 prod_{m>1}(phi1-phi_m)eps^3 / prod_j(1 - kappa z_j) = O(1).
  The task's first-order formula n'D_kappa n'^{-1} ~ eps[diag(c - phi) + (phi_i - phi_j) n'_ij] is wrong for the (1,4) entry:
  the n'_12, n'_13 paths cancel it (sum_j r_j phi_j = 0).  => a common rotation on <0| would also give only x (no xy):
  the singlets need p23 AND p32, impossible from <0| with any rotation that is a polynomial (deg <= 2) in the twist on A.
* 2603-type twist-DRESSED references do reach generic polarisations: dress the pair-1 A site by (zeta - Gc^T)^{-1} and the F site
  by (zeta' - Gc)^{-1} (diag frame: v_j ~ lambda_j/(zeta - z_j), u_j ~ 1/(zeta' - z_j) -- Cauchy-like); with zeta = 1 + eps xi,
  zeta' = 1 + eps xi' the physical images are O(1) with ALL components (sympy, eps -> 0):
     v ~ (prod_{m!=1}(xi-phi_m)/prod_{m!=1}(phi_1-phi_m)-type, ..., 1),   u ~ (1/prod_m(xi'-phi_m)-type, ..., 1)
  (exact expressions in t17_frames.py output), so p23 != 0 != p32.  BPS (v.u = 0) iff prod_j(zeta - z_j) = prod_j(zeta' - z_j)
  (v.u = [1/((zeta-x)(zeta'-x))][z1..z4] = (1/c(zeta) - 1/c(zeta'))/(zeta' - zeta), c(x) = prod(x - z_j)): a one-parameter
  family of BPS pair-1 polarisations seeing the singlets.  These dressings are NOT in the localised word span (they need companion
  A components e2, e3 and F component e1), so they must enter as the reference (2603's step (i): prove <W|Psi> is a Q-bracket
  determinant).  SoV word expansions of such dressed references are dense (t15_dressed.py: 11-25 words for one dressed site, 76-256
  for several) -- as expected, the determinant structure lives in the measure, not in sparsity.

## 6. SoV norm exactly at the ABJM point: complete regularised bra AND ket bases (L = 2, exact mod p)
* Standard bra and ket word bases both have rank 81 at thA = thF + 2h (t11_ketrank.py).
* Exact ket identities (t12_ketident.py): t2(thA) = t1(thA) t1(u0), t3(thA + h) = 0 (mirror of t2(u0) = t1(u0)t1(thF), t3(u0) = 0).
* Regularised ket basis t1(th1)^m1 K2 t1(th3)^m3 K4 |0bar>, K in {1, t1(u0), Wt = t2'(thA) - t1(thA) t1'(u0), Zt = t3'(thA+h)}:
  rank 256 (t14_ketvariants.py).  Other splits: t2' - t1'(thA)t1(u0): 241 (240 with a k=3 exception); full derivative
  t2' - t1'(thA)t1(u0) - t1(thA)t1'(u0): 248 (256 with the k = 3 exception t2(thA) for m = 3).
* With the ABJM-point agent's regularised bra basis (rank 256 confirmed here) both bases are complete, so
  1/N = sum_{x,y} tau~_y M_yx tau_x (M = Gram^{-1}) holds identically at the ABJM point; t18_norm_abjm.py records the Gram structure
  and checks it on eigenstates: t18 result: ranks 256/256; Gram 3323 nonzeros, charge sums 12..20 (generic-theta rule
  S_x + S_y = 12 and triangularity LOST at the ABJM point); M = Gram^{-1} 4267 nonzeros, cond 4.2e10;
  1/N = sum tau~ M tau vs direct: 1.4e-5 (float64, 29 states); tau_x = prod_pairs lambda(a,k) tauF^k: 2e-13..1e-4.

## Status
- [x] Task 1 at L = 2 exact (negative: K_kappa not a product on <0|); 16-word regularised scan running (t02)
- [x] Task 2: physical identification (negative for the even bridge; dressed references identified as the way forward)
- [x] Task 3: formula + untwisting verification L = 2 (mp); L = 3 untwisting NOT done (mp cost)
- [x] Task 4: Q-form (tauF, x, w, y Cramer formulas verified); norm via regularised bases
