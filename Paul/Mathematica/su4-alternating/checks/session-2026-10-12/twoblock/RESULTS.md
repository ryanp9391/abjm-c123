# Two-block rotations of the companion SoV vacuum, L = 4 (A F A F) — results (2026-10-12)

Notebook: `Paul/Mathematica/su4-alternating/sov-twoblock-expansion-L4.wb` (Wolfram cells UNRUN — kernel died / licence
expired; python results pasted as markdown under each code cell). All numbers below are python: exact = mod-p over >= 2 primes
(p ~ 2^26) or CRT + rational reconstruction (verified on an extra prime); sympy symbolic; mpmath 60/90 digits; float = numpy.
Point: h = 1/2, theta = (2/7, -3/5, 8/11, -11/8), z = (3/2, 5/4, 2/3, 4/7).

## Setup
- t1, t2, t3 (quantum minors) commute [exact]. c1 = 603680/25947. Calibration (i) <Om|t1(th1)^-1 = c1 <Om|rho1(G^-T) holds for
  GLOBAL BPS rotations only (two-block BPS: fails, 147 nonzero residual comps) [exact]. (ii) twisted-vacuum eigenbra [exact].
- Word basis rank 256, exactly z-independent [exact]. <x|Psi> = tau-product <0|Psi> on all 256 twisted states [float 3e-9].

## SoV-charge grading [exact]
- M (words x product basis) is block diagonal in |x| = a1+k2+a3+k4 = sum(A: i-1, F: 4-j); blocks 1,4,10,20,31,40,44,40,31,20,10,4,1;
  nnz(M) = nnz(M^-1) = 2929. Dual basis Y = M^-1 = coefficient tensor C_x^{ijkl}; c_x(B) = (B.Y)_x.
- Operator-Schmidt ranks of Y: {12}|{34} 169/256, site 1 13, site 2 16, site 3 16, site 4 13, {13}|{24} 201.

## Two-block expansions [exact] (# words, 16x16 block rank)
- generic two-block / global BPS: 256, 16. block{1,2} generic + block{3,4} = SoV vac: 36 (a3 = 0, k4 <= 2), 3.
  block{3,4} generic: 130, 16. block{4,1}: 96, 13. block{2,3}: 83, 10. pair-1 common g = 1-kG^T: 20, 3; pair 2: 52, 9.
- flags e1 p(G^T), e4 q(G), deg p + deg q <= 2 (= Lagrange-weighted twist-adapted lam_j p(z_j), q(z_j)): 1..10 words on block 1.

## Parameter fixing (main result) [sympy exact + 60/90 digits]
- <0| (diag frame) = (lam x 1)^{x2}, lam_j = 1/prod_{m!=j}(z_j - z_m); <0| = lam1^2 (vac x vac) rho(n'D'), r = lam/lam1,
  n' = ((1,1,1,r4),(0,1,0,-r2),(0,0,1,-r3),(0,0,0,1)), D' = diag(1,r2,r3,1): finite as z -> 1.
- Block rotations giving the physical bra: g1 = 1 + b E, g2 = 1 + b' E (mod stabilisers, common Borel), E = E14 in the twist
  eigenbasis (not a twist polynomial): <0|rho12(g1)rho34(g2) = lam1^2 (W_b x W_b') rho(n'D'); for HW Psi0 of weight (w1,0,0,-w1):
  overlap = lam1^2 (b'-b)^{2-w1} C0, C0 = <W_-1 x vac|Psi0> = <vac x W_-1|Psi0>.
- Sparsest: b' = 0: B(b) = (e1 + b lam1 v_z4) x (e4 - b l1) x e1 x e4 (companion), v_z4 = e1 prod_{m!=4}(G^T - z_m),
  l1 = lam1 e4 prod_{m!=1}(G - z_m). 36 SoV words, a3 = 0, k4 <= 2; b^0/b^1/b^2 pieces 1/17/36 words (vacuum/adjoint/singlet).
- Untwisting (theta = theta_ABJM + eps theta0, z = 1 + eps phi, phi = (.9,.4,-.3,-1.1)): C_pred = <B(b)|Psi_eps>/(lam1^2 (-b)^{2-w1})
  -> data C2 linearly in eps, b-independent: singlets -0.15359/-0.15766/-0.15807 vs -1/sqrt40, 0.19981/0.20369/0.20408 vs 1/sqrt24
  (eps = 1e-2/1e-3/1e-4); vacuum 1; [1,0,1], [0,2,0] -> 0.
- 36-term SoV sum = direct overlap to 1e-44 (60 dig, eps = 1e-2) and 1e-65 (90 dig, eps = 1e-3); <0|Psi>/lam1^2 ~ eps^2 ([0,2,0]),
  eps^3 (adjoints), eps^6 (singlets). Float64 fails (cond M 1e18-1e25).

## Determinant structure
- K_b = sum_a tau_a(th1-h) P_a(tau1(th2), tau1(th4)) [float 7e-10 on 256 states] = det[Q_i(th1 - beta h) (i=1..3, beta=0..3);
  (-1)^{beta+1} P_beta] / Q123(th1-h) (Laplace, given dictionary tau_a = D_a/Q123).
- Pair matrices ranks 4,2,1 (b^1), 4,4,3 (b^2); site ranks 4,4,1,3; P_a non-factorisable; CP rank 7 (ALS, generic).
- Vandermonde ansatz (coordinator), two-site bilinear forms K = g^T Kmat g' [exact]:
  * F-F (sites 2,4; u = gamma e4 + delta e3): K = (g2+d2 rho2)(g4+d4 rho4) - (h/th24) g2 d4 rho2 - (h/(th24+h)) d2 d4 rho2 rho4,
    rho2 = tau1(th2)/(h th12 th23 (th24-h)), rho4 = tau1(th4)/(h th14 th24 th34): EXACT as written (rho2 in the g2 d4 term).
    The reading mu = (x2-x4)/(x2-th4) (needs rho4 there) is FALSE. det Kmat = h rho2 [rho2/th24 - rho4/(th24+h)] not separable
    => no Vandermonde x separable reading and no determinant with one bra-carrying row per site (any constant rows) exists.
    u in span(e2,e3,e4): 3x3 form generic rank 2.
  * A-A (site 1 {e1,e4}, site 3 {e1,e2}): det != 0, not separable => excluded. README {1,3} formula: det = (h/th13)(1 - o1 o3), excluded.
  * A-F (A1,F2),(A1,F4),(A3,F2),(A3,F4) in these directions: det = 0 exactly => K factorises, no A-F cross factor.
    But single-site factors are nonlocal (site 3 + e2 involves tau1(th2), tau1(th4), tau1(th1-h); site 4 + e3 involves tau1(th2)),
    and site 1 + e2 gives K ~ 43590625 r2 - 4877568 r4 - 72678508 tau1(th1-h): A-F cross terms in other directions.
  * Physical bra uses all directions -> no block determinant of the conjectured form in the left-basis/tau language.
  * Diagnosis: Y (B-eigen right basis) gives unit measure and non-separable right wave functions; a 2603/GLMR determinant needs a
    C-eigen right basis with non-diagonal measure. Commuting candidate (reflected RV B, reversed minors) found numerically but its
    eigenvalue roots are not on theta + hZ (not SoV-type). Open.

## Scripts (this folder)
core.py (chain, t_a, words; float & mod-p), exact.py (CRT/ratrec), wxf.py (WXF reader), states.py, tbtools.py, nb_common.py,
nb_cell*.py (notebook result ports), t01-t30*.py (exploration), s01-s05*.py (Vandermonde ansatz tests). Run from this folder;
some need cached .pkl files produced by t05/t08/nb_cell15 (dual_mod*.pkl, Y_mod20.pkl).
