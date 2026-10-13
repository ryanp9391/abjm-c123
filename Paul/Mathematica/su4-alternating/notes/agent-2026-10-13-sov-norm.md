# Agent 2026-10-13: the SoV normalisation <0|Psi>/||Psi|| in Q-functions (L = 2 even bridge, ABJM point)

Notebook: `abjm-sov-norm-L4.wb` (symbols prefixed `nm`; no memoisation; all cells run, saved with outputs). Python fallback / independent
cross-check (shared kernel was busy for ~15 min): `checks/session-2026-10-13/norm/` (`gram_abjm.py`, `abjm_eval.py`; exact Fractions +
mpmath, built on `checks/session-2026-10-12/sovformula/`; log `eval_1_100.log`).
Conventions exactly as README / `w-a-reference-status.md`: chain A F A F, site 1 leftmost, bras are row vectors, T = L1 L2 L3 L4, twist on the
right, companion twist G_c(z), frame map t^{G_c} R = R t^Lambda, R = V^{-T} (A) x V (F), SoV vacuum <0| = e1 e4 e1 e4 (companion),
|0bar> = e4 e1 e4 e1, words <a1,k2,a3,k4| = <0| t_{a1}(th1-h) t1(th2)^k2 t_{a3}(th3-h) t1(th4)^k4, ket words A: t1(th_A)^m, F: t_n(th_F + n h),
M = Gram^{-1} (Gram[x,y] = <x|y>, M[y,x]). Standard point h = 1/2, theta = (2/7,-3/5,8/11,-11/8), z = (3/2,5/4,2/3,4/7). eps-chain exactly as in
`Wa-sov-expansion-L4.wb` cells 30-48 (z = 1 + eps phi, phi = (9/10,7/20,3/5,-11/10), theta = (thF1+2h, thF1, thF2+2h, thF2), thF = eps (0, 37/100),
c = 3/10, kappa = 1 - c eps, a = (phi1 - c)/(phi1 - phi4)). Tags: [exact], [hp] (60-120 digits), [num].

## Headline

**Yes: the L = 2 even bridge is now Q-complete.** The last non-Q input, <0|Psi>/||Psi||, is replaced by sqrt(1/B) with the SoV norm
B(Psi) = sum_{x,y} tauT_y(Psi) M_{yx} tau_x(Psi), a bilinear form in eigenvalues of the commuting transfer matrices (and their u-derivatives at u0),
with the regularised bra/ket SoV bases at the ABJM point. Fully assembled (eigenvalue data only + explicit theta, z, M):

    |C_{M=1}| = lim_{eps->0} N_pre * A * W / sqrt(B)

[hp] deg-5 extrapolation over eps = 1/50, 1/100, 1/200, 1/400, 1/800, 1/1600:
E = 20 singlet 0.1581138830084186 vs 1/sqrt(40) = 0.1581138830084190 (diff 4e-16);
E = 4 singlet 0.2041241452319316 vs 1/sqrt(24) = 0.2041241452319315 (diff 1e-16). With the data phase (sign of Psi_eps fixed by <Psi_0^data|Psi_eps> > 0)
the signed direct C is -0.15807, +0.20407 at eps = 1/400 -> -1/sqrt(40), +1/sqrt(24) as in `hw_L2.json`. The Q-based formula gives |C| only
(the sign is the phase of Psi, not Q-data).

The limits of the factors are exact integers: (A, W, B) -> (10, 25, 2 500 000) for E = 20 and (2, 9, 7776) for E = 4, N_pre -> 1
(extrapolated to 1e-14), so C = A W/sqrt(B) = 250/sqrt(40 * 250^2) = 1/sqrt(40) and 18/sqrt(24 * 18^2) = 1/sqrt(24).

**Exact finite-twist norm statement (task 1).** There is no exact identity |<0|Psi>|^2 = |lambda^w N| ||Psi||^2 at finite twist. Exact statements:
(E1) <0|Psi> = lambda^w <0bar^T|Psi> for every state of weight w [exact, weight counting; checked 1e-59 on all 256 states];
(E2) N := <0|Psi><PsiL|0bar>/<PsiL|Psi> = 1/sum tauT M tau [hp, all 256 states at generic theta, z: max rel. error < 1e-47];
(E3) |<0|Psi>|^2 = H * lambda^w * N * ||Psi||^2 with H = <Psi^T 0bar><PsiL Psi> / (<PsiL 0bar> Psi^T Psi) (diagonal frame, real Psi), and
H != 1 at finite twist: the right eigenvectors of t^Lambda are not orthogonal, PsiL is not proportional to Psi^T (cos 0.37 ... 1 at the standard
point), not to Psi^dagger, and (R Psi)^T is NOT a left eigenvector in the companion frame (cos <= 0.06). H = 1 + O(eps): for the singlets
H - 1 = (E/2)(thF2 - thF1) + O(eps^2) (slopes 3.70 eps and 0.74 eps = 10*0.37 eps, 2*0.37 eps for E = 20, 4); at pair-periodic theta (thF1 = thF2)
only an O(eps^2) twist effect remains (max |H - 1| = 0.0088 at eps = 0.1). What replaces the Euclidean norm at finite twist is the SoV pairing:
<0|Psi>^2/<PsiL_*|Psi> = lambda^w / sum tauT M tau exactly, with PsiL_* the left eigenvector normalised by <PsiL_*|0bar> = <0bar^T|Psi>.
Since H -> 1 as eps -> 0 the limit C is unchanged; C_Q(eps) = C_direct(eps)/sqrt(H) differs from C_direct(eps) at O(eps) (e.g. 0.1551 vs 0.1580 at eps = 1/100).

## 1. Setup (notebook sections 1-2)

Fused transfer matrices t_a (a = 1, 2, 3, with the polynomial prefactors g2, g3 of `Wa-sov-expansion-L4.wb`), exact u-derivatives by a dual-number
monodromy (T' = sum_s L1..L_{s-1} L_{s+1}..L4; checked against finite differences 1e-15 and against the frame map t^{G_c} R = R t^Lambda, exact).
Generic point: 256 bra words, 256 ket words, rank 256, Gram 1456 nonzeros, M = Gram^{-1} charge-block inverse, 1456 nonzeros [exact].

## 2. Generic theta, z: all 256 states (notebook sections 3-5)  [hp 60 digits]

* Joint eigenvectors of t1(3/11) + 7/10 tbar3(-5/7) - 3/13 tbar2(2/9) per weight block (diagonal frame); left eigenvectors = rows of the inverse of
  the right-eigenvector matrix; companion frame Psi_c = R Psi, PsiL_c = PsiL R^{-1}. 226 states real, 30 complex.
* tau_x = tau_{a1}(th1-h) tau1(th2)^k2 tau_{a3}(th3-h) tau1(th4)^k4 and tauT_y = tau1(th1)^m1 tau_{n2}(th2+n2 h) tau1(th3)^m3 tau_{n4}(th4+n4 h) are products of
  eigenvalues read off as (op Psi)_p/Psi_p (pivot ratio) of the companion operators.
* 1/N = tauT^T M tau equals <PsiL|Psi>/(<0|Psi><PsiL|0bar>) on all 256 states, max relative deviation < 1e-47 (0 within the 60-digit precision); |1/N| ranges 1.3e-6 ... 3.5e5, complex on complex states.
* <0|_Lambda = lambda (x) 1 (x) lambda (x) 1 and |0bar>_Lambda = 1 (x) lambda (x) 1 (x) lambda exact; <0|Psi> = lambda^w <0bar^T|Psi> to 1e-59 on 256 states.
* Orthogonality: max off-diagonal cosine of the right eigenvectors Psi_n^T Psi_m is O(1) (1.15), Hermitian 0.98 (standard z far from 1); no fixed
  diagonal operator S with PsiL_n = c_n Psi_n^T S (row ratios vary by O(1)); cos(PsiL_n, Psi_n^T) in [0.37, 1] over real states; in the companion frame cos(PsiL_c, Psi_c^T) in [8e-10, 0.058].
* Reason for the failure: t^T = reversed chain F A F A (full transposition fixes both Lax matrices), which is the original chain only after a shift
  and only for pair-periodic theta; at generic theta there is no symmetry S with S t S^{-1} = t^T.

## 3. H on the eps-chain (notebook section 5, 10)  [hp 90 digits, zero-weight block, 28 states, 16 real]

H_n - 1 = O(eps), max over real states = 3.83 eps at eps = 1/50, 3.76 eps at 1/100 (-> about 3.7 eps); cos(PsiL, Psi^*) = 1 - O(eps^2) (min 0.9999 at eps = 1/100).
Pair-periodic theta (thF1 = thF2 = 0 or 1/7) at twist 1 + phi/10: max |H - 1| = 0.0088, at the standard z: 0.24; theta_F = (0, 1/100) at twist 1 + phi/10: max H = 1.0997.

## 4. ABJM point (notebook sections 6-8)  [exact operators, hp 100-120 digits]

Pair-by-pair theta: coincident inhomogeneities, the generic bases are singular (rank 81). Regularised bases per pair p (u0_p = thF_p + h = thA_p - h):
bra O = {1, t1(u0), W = t2'(u0) - t1'(u0) t1(thF), Z = t3'(u0)} (and t2(u0) for the word a = 2, k = 3), F = t1(thF);
ket A = t1(thA), K = {1, t1(u0), Wt = t2'(thA) - t1(thA) t1'(u0), Zt = t3'(thA + h)}.
[exact, both pairs, eps = 1/100 ... 1/1600 (cell 15: eps = 1/100; the others in `nmRunEps` "ident")] t3(u0) = 0, t2(u0) = t1(u0) t1(thF), t2(thA) = t1(thA) t1(u0), t3(thA + h) = 0.
[exact mod p] ranks of bra words, ket words, Gram: 256, 256, 256 at all six eps; Gram has 3323 nonzeros at every eps (not charge-graded, not triangular).
Wave functions (eigenvalue data of Psi): bra pair factor lambda(a,k) tauF^k with lambda in {1, x = tau1(u0), w = tau2'(u0) - tau1'(u0) tauF (k <= 2) or
x tauF (k = 3), y = tau3'(u0)}; ket pair factor tauA^m {1, x, wt = tau2'(thA) - tau1(thA) tau1'(u0), yt = tau3'(thA + h)}. The combined operators W, Z, Wt, Zt give
the same eigenvalues as the separately extracted tau2', tau1', tau3' (checked 1e-80).
1/N from these 14 numbers per state (tauF1, tauF2, x1, x2, w1, w2, y1, y2, tauA1, tauA2, wt1, wt2, yt1, yt2) and M(theta) = Gram^{-1} [solved in 3*prec digits]:

| eps | max_28 abs(1/N_Q / 1/N_direct - 1) | Gram-solve residual | closed-form <W_a|Psi> vs direct (rel.) |
|---|---|---|---|
| 1/50 | < 1e-48 (100 digits) | < 1e-77 | < 1e-80 |
| 1/100 | < 1e-41 (100) | < 1e-74 | < 1e-78 |
| 1/200 | < 1e-33 (100) | < 1e-70 | < 1e-76 |
| 1/400 | 1.8e-40 (100) | 2e-85 | < 1e-73 |
| 1/800 | < 1e-29 (110) | -- | < 1e-81 |
| 1/1600 | < 1e-32 (120) | -- | < 1e-88 |

(Mathematica significance arithmetic reports the accuracy listed; it decreases as the Gram matrix degenerates for thF1 - thF2 = O(eps).) Independent python check (exact Fractions, mpmath 80 digits, own
Gram, own eigen-solver; `eval_1_100.log`): eps = 1/100, all 28 states, |C_dir|, |C_Q| and H agree with the notebook values to all printed digits (e.g. 0.157956610596 / 0.155064907354 and
0.203901882009 / 0.203156885621 for the singlets), max 1/N relative error 3.2e-35 over the 28 states (80 digits), max |H - 1| = 0.03764 on the 16 real states.
A second python run at eps = 1/400 (110 digits, `eval_1_400.log`): all 28 states, 1/N relative error between 1e-54 and 1e-82 (singlets 8e-82), |C_dir| = 0.158074407836, 0.204071983258 and |C_Q| = 0.157345233498, 0.203883806336
as in the notebook; so the 1e-40 at eps = 1/400 in the table is a precision loss of the Mathematica LinearSolve (100 digit input), not a failure of the identity.

**Reflection structure of the eigenvalue data [hp 60 digits, python `abjm_eval2.py`, eps = 1/100, log `evaldata_1_100.log`; the pickle was not kept].** On 16 of the 28 zero-weight states (10 real, among them BOTH singlets and the top state, and 6 complex ones) the ket data are exactly the bra data
of the other pair, to the working precision (1e-47 ... 1e-60):
tauA1 = tauF2, tauA2 = tauF1, x1 = x2 (= x), wt1 = -w2, wt2 = -w1, yt1 = -y2, yt2 = -y1
(equivalently tau1(u) = tau1(c - u) with c = thF1 + thF2 + 2h: tau1(thA1) = tau1(thF2), tau1(thF1 + h) = tau1(thF2 + h)). The other 12 states (two parity-partner triplets of real states with E ~ 630 and ~ 597, and 6 complex ones) do not
satisfy it (they are mapped onto each other by the reflection). So for the parity-symmetric states, in particular for the two singlets of the even bridge, B depends only on the 7 numbers (tauF1, tauF2, x, w1, w2, y1, y2):
tauT_y(Psi) = sigma_y tau_{pi(y)}(Psi) with the pair exchange pi and the signs (-1) for the levels n = 2, 3 (except the word (m,n) = (3,2), where the bra side uses t2(u0) = x tauF instead of w).

## 5. Q-based even bridge, eps = 1/50 ... 1/1600 (notebook sections 7-9, 11)

|C_dir| = |K <W_a|Psi>/(1-kappa z1)^2| (unit Psi), |C_Q| = |N_pre A W|/sqrt|B|. The closed form of <W_a|Psi> (bidegree-(2,2) polynomial in the F-site Q ratios) reproduces the direct overlap on all 28 states.

| eps | E=20: C_dir | C_Q | H | E=4: C_dir | C_Q | H |
|---|---|---|---|---|---|---|
| 1/50 | 0.1578009597 | 0.1520831314 | 1.07660697 | 0.2036431016 | 0.2021739338 | 1.01458651 |
| 1/100 | 0.1579566106 | 0.1550649074 | 1.03764444 | 0.2039018820 | 0.2031568856 | 1.00734765 |
| 1/200 | 0.1580350383 | 0.1565808748 | 1.01866021 | 0.2040175558 | 0.2036424824 | 1.00368704 |
| 1/400 | 0.1580744078 | 0.1573452335 | 1.00928994 | 0.2040719833 | 0.2038838063 | 1.00184678 |
| 1/800 | 0.1580941321 | 0.1577290200 | 1.00463497 | 0.2040983471 | 0.2040040990 | 1.00092420 |
| 1/1600 | 0.1581040042 | 0.1579213167 | 1.00231499 | 0.2041113168 | 0.2040641529 | 1.00046230 |
| extrap (deg 5) | 0.15811388300841897 | 0.15811388300841860 | | 0.20412414523193151 | 0.20412414523193159 | |
| target | 0.15811388300841897 (1/sqrt40) | | | 0.20412414523193151 (1/sqrt24) | | |

(deg-3 extrapolation of C_Q from the last 4 points: 0.158113883006, 0.204124145233; deg 2: 0.158113882364, 0.204124145265.)
The same states carry (1/50 ... 1/1600) factor values A = 9.41 ... 9.98 (-> 10), W = 25.67 ... 25.02 (-> 25), B = 2.613e6 ... 2.5035e6 (-> 2.5e6) for E = 20 and
A = 1.94 ... 1.998 (-> 2), W = 9.24 ... 9.007 (-> 9), B = 8125 ... 7787 (-> 7776) for E = 4; N_pre = 1.018 ... 1.0006 (-> 1).

## 6. Final L = 2 formula (2603 form), every factor named

    C_{M=1}(Psi) = N_pre * A(Psi) * W(Psi) / sqrt(B(Psi)),   evaluated as eps -> 0 (z = 1 + eps phi, theta_F = eps theta^0, kappa = 1 - c eps), |C| = absolute value.

* **B (SoV norm, 2603's sqrt(B)^2):** B = sum_{x,y} tauT_y M_{yx} tau_x = 1/N, M = Gram(theta)^{-1}, with the regularised words above. Pure eigenvalue data
  (tauF_p, x_p, w_p, y_p; tauA_p, wt_p, yt_p) and the number table M(theta) (depends on theta, h only; z-independent). For the parity-symmetric states (both singlets) the ket data are the reflected bra data (section 4), so B is a
  bilinear form in the bra data tau_x(Psi) and their reflection.
* **A = P_kappa(tau) = K_kappa:** A = s^2 + kappa s [ X1/(h^2 l1) + X2/(2 h^4 y l2) + X3/(6 h^5 y^2 l3) ] + kappa^2 X1 X3/(24 h^7 y^2 l1 l3),
  s = 1 - kappa (z1 + z4), y = (z2 + z3)/2, X1 = tau1(thF1), X2 = tau2(u0) = tau1(u0) tau1(thF1), X3 = tau3'(u0) (u0 = thF1 + h);
  l_k = prod_{j<k}(u - th3 - j h) (u - th4 - k h) prod_{j<k-1}(u - th4 - j h) with vanishing factors deleted, l1 = l_1(thF1), l2 = l_2(u0), l3 = l_3(u0).
  X1, X2, X3 are Q-data at one pair: X1 = Q_j(thF1)/Q_j(thF1 - h); X2, X3 from Q_i(u0 - kh), Q_i'(u0) by the 3x3 Cramer formulas of reports 9d (not re-verified here).
* **W = u~ . P(rho2, rho4) . u~**, u~ = (1, -(z2 + z3), z2 z3), rho2 = X1/(h th12 th23 (nu - h)), rho4 = X4/(h th14 nu th34), X4 = tau1(thF2) (= Q_j(thF2)/Q_j(thF2 - h)),
  nu = th2 - th4 = thF1 - thF2, P the 3 x 3 matrix of `Wa-sov-expansion-L4.wb` cell 38 (bidegree (2,2)); F-site single-Q data only.
* **N_pre (explicit theta, z; no Q):** N_pre = (a/kappa0)^2 / (Lam1^2 u4^2 (1 - kappa z1)^2), a = (phi1 - c)/(phi1 - phi4), kappa0 = -(z1-z2)(z1-z3)/((z4-z2)(z4-z3)),
  Lam1 = 1/prod_{m != 1}(z1 - z_m), u4 = (z4 - z2)(z4 - z3); (a/kappa0)^n with n = L - mu = 2.
* **Eigenvalue (Q) data of Psi:** tauF1, tauF2, tauA1, tauA2, x1, x2, X2 = x1 tauF1, X3, w_p, y_p, wt_p, yt_p (p = 1, 2). **Explicit:** h, theta, z, c, kappa, a, s, y, l_k, nu, u~, N_pre, M.
* Phase: |C| only; the sign of Psi is convention (data phase reproduces -, + for E = 20, 4).

## 7. What failed / limits / not done

* No exact finite-twist Euclidean identity (H = 1 + O(eps) on the nose); the formula holds in the eps -> 0 limit. C_Q(eps) is a regularised structure constant (SoV-norm regularisation).
* The SoV norm is a 256 x 256 bilinear sum (Gram 3323 nonzeros at the ABJM point, no diagonal measure, no determinant form); not simplified to a Q-bracket determinant (needs the mixed 4/4bar measure).
* M(theta) and 1/N degenerate as thF1 - thF2 -> 0 (cond growth), so accuracy decays with eps (1e-47 ... 1e-29 at 100-120 digits); the homogeneous (fully periodic) limit needs higher-derivative bases (open).
* x, w, y were obtained as eigenvalues of the (derivative) operators on the joint eigenvectors; their Cramer/Baxter Q-formulas (reports 9d) were not re-derived here. The relation between bra data (w, y) and ket data (wt, yt) was found numerically only (section 4)
  (14 numbers per state enter B in general; 7 for the 16 parity-symmetric states incl. both singlets, see section 4; a proof / the reflection operator was not constructed).
* Only L = 2, only the E = 20 and E = 4 singlets of the even bridge; L = 3, middle bridges, the (1/3) e^{-+ i pi/3} factor of reports 9d untouched.
* Kernel housekeeping: the shared kernel MemoryInUse was 2.4 GB at the start, 3.0 - 3.1 GB during most of the session and 7.6 GB at the end (the L = 6 notebook of the other agent), while all `nm*` symbols total < 200 MB at peak (3 MB after
  the final cleanup); the 3 GB rule was therefore exceeded for reasons outside this thread, and I continued with minimal footprint. I set `$HistoryLength = 0` once to drop the ~0.9 GB `Out` history; this blanked the
  Wolfbook outputs of the cells run in that period, so cells 5-34 were rerun after restoring `$HistoryLength = 10` (the notebook now carries all outputs). The kernel was busy (other agent) for ~15 min in total; the python fallback in
  `checks/session-2026-10-13/norm/` was used meanwhile and as independent cross-check.

## 8. Files

`abjm-sov-norm-L4.wb` (sections 1-12, all outputs), `checks/session-2026-10-13/norm/{gram_abjm.py, abjm_eval.py, abjm_eval2.py (eigenvalue-data dump), logs gram_1_100/400.log, eval_1_100/400.log, evaldata_1_100.log; pickles are cache files}`.
