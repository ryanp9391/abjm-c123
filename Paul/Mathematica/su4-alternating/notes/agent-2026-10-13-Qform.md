# Agent 2026-10-13: consolidated Q-form notebook and the eps -> 0 limit of A and W at L = 2

Notebook: `even-bridge-Qform.wb` (prefix `qf`, no memoisation, all cells run, saved with outputs). Conventions exactly as README / `w-a-reference-status.md`.
Scope (after the coordinator's change): the SoV norm B is the known Q-function norm (Paul) and is **not** rederived or limited here; only
A = P_kappa and W = B_{z2} B_{z3} are taken to the untwisted limit. Tags: [exact] rational/algebraic, [hp] fixed high precision, [num] double.

## What the notebook contains

* Part 0: one set of builders (Module-local sparse monodromy with exact u-derivative at L = 2; vector propagation at L = 3), checked: frame map
  t^{G_c} R = R t^Lambda (k = 1,2,3), <0|_Lambda = lambda x 1 x lambda x 1, derivative vs central difference, vector propagation vs matrices [exact].
* Part 1 (a) <W_a^L|K_kappa = (1 - kappa z1)^2 <W_a' x W_a^{L-1}|, six (a, kappa), z2 = z3, L = 2 [exact]; control z2 != z3 fails (rel. residual 0.069).
  (b) <0_w|Psi> = (v1 u4)^L <W_a^L|Psi>, w = (w1,z2,z3,w4): all five (mu,0,0,-mu) sectors, two w1, 40 digits [hp]; off-sector controls; Moebius map symbolic [exact].
  (c) <W_a^L|Psi>/<0|Psi> = (a/kappa0)^n (lambda1 u4)^-L B_{z2}(y) B_{z3}(y), y_s = tau_1(theta_s)/kappa_s: L = 2, all 108 states of sectors n = 0..4 at generic and ABJM-pair theta, < 1e-41 [hp 50 digits]; L = 3, sectors n = 0,1,2 (1+3+3 states), < 1e-49 [hp 50 digits].
* Part 2 (a): |C_{M=1}| = lim N_pre A W |<0|Psi>|/||Psi|| on the L = 2 eps-chain (eps = 1/100..1/1600, 90-110 digits): Q-form of <W|Psi> equals the direct overlap to < 1e-70 on all 28 zero-weight states; |C| -> 1/sqrt40, 1/sqrt24 (degree-4 extrapolation, ~1e-15); (A, W, N_pre) -> (10, 25, 1), (2, 9, 1); |<0|Psi>|/||Psi|| -> 1/sqrt(2.5e6), 1/sqrt(7776).
  (b) kept from the first pass, not extended: regularised-basis B = sum tauT M tau (eps = 1/50..1/1600), B -> 2.5e6, 7776. (c) L = 3 table with the norm from the eigenvector.

  (c) result [hp 70 digits, eps = 1/100..1/800]: Q-form (norm from the eigenvector) = direct to ~1e-41; extrapolated |C| = 0.17112641, 0.027137747, 0.045900846 (x2), 0.15348408 (x2) = data. At L = 3 the factors are NOT separately finite: W ~ eps^2, |<0|Psi>|/||Psi|| ~ eps^-3, N_pre ~ eps (only the product is); A -> 3.65685, 7.65685, 14.658, 1.842 (the known K_1 values).

## Part 3 results (L = 2, ABJM point z = 1, theta_F = (0,0), theta_A = (1,1), h = 1/2)

**3.1 [exact]** Untwisted singlets (2-dim joint kernel of the six raising/lowering generators in the zero-weight sector). Per pair:
E = 20: tau_F = tau_1(0) = 1/4, tau_1'(0) = 1, x = tau_1(1/2) = 3/4, tau_A = 1/4, w = 3/4; E = 4: tau_F = -1/4, tau_1'(0) = 1, x = 1/4, tau_A = -1/4, w = 1/4.
The *plain* derivatives y = tau_3'(u0) and y~ = tau_3'(3/2) vanish exactly: at theta_F1 = theta_F2 the prefactor of t_3 has a double zero at u0 (both F sites), so
t_3'(u0) = 0 and the regularised hat t_3(u0) of K_kappa is t_3''(u0)/2. At finite nu = theta_F1 - theta_F2, y_p = O(nu), divided by the localisation constant lambda_3 ~ nu.

**3.2  A -> K_1 [exact + hp].** With A = tau_F/(h^2 lambda_1), B = tau_2(u0)/(2 h^4 y lambda_2), C = hat tau_3(u0)/(6 h^5 y^2 lambda_3):
P_kappa = s^2 + kappa s (A + B + C) + kappa^2 A C/4 -> K_1 = 2 - A - B - C at kappa = 1, y = 1, s = -1 (lambda -> (1/2, -1/8, 3/8) with vanishing factors dropped).
Exact singlets: (A, B, C) = (2, -12, 2) and (-2, 4, -2), hat tau_3 = -+ 9/64, A C = 4, K_1 = 10 (E = 20) and 2 (E = 4). The normalised A, B, C of the twisted eps-chain extrapolate to these values (hp, ~1e-12).
With tau_2(u0) = x tau_F and A C = 4 (translation-invariant states, X3^ X1 = (9/256) P):
**K_1 = 2 - 8 tau - 1/(2 tau) + 64 x tau,  tau = tau_1(theta_F), x = tau_1(theta_F + h)** [exact on the singlets].

**3.3  W -> B_1^2 [analytic derivation; exact on singlets; hp on all 28 states].**
kappa_2 = 2 h^2 nu (nu - 2h), kappa_4 = 2 h^2 nu (nu + 2h), g = nu^2/(nu^2 - h^2).
* The normalised y_s have **no finite untwisted limit**: nu y_2 -> -tau/(4 h^3), nu y_4 -> +tau/(4 h^3) (y_{2,4} ~ -+2 tau/nu at h = 1/2).
* Finite combinations: g y_2 y_4 -> tau^2/(16 h^8); y_2 + y_4 -> -(tau + h tau')/(4 h^4), with tau' = d tau_1/du at theta_F (the poles cancel because tau_1(theta_F1) - tau_1(theta_F2) = nu tau' + O(nu^2); the O(eps) twist shifts cancel in the difference).
* Hence B_{z2}, B_{z3} -> B_1 = 1 + (tau + h tau')/(4 h^4) + tau^2/(16 h^8) = **1 + 4 tau + 2 tau' + 16 tau^2** (h = 1/2) and W -> B_1^2.
* Singlets [exact data]: E = 20: tau = 1/4, tau' = 1 -> B_1 = 5; E = 4: tau = -1/4, tau' = 1 -> B_1 = 3. The eps-chain numbers (nu y_2, nu y_4, y_2+y_4, g y_2 y_4, B_1, B_z2, B_z3, sqrt|W|) extrapolate to (-1/2, 1/2, -3, 1, 5, 5, 5, 5) and (1/2, -1/2, -1, 1, 3, 3, 3, 3) [hp].
* Parity shortcut: if tau_1(u) = tau_1(1 - u) is a degree-4 polynomial with leading coefficient 4, tau'(0) = 4x - 4 tau - 1 and B_1 = 16 tau^2 - 4 tau - 1 + 8 x. Exact on the two singlets (5, 3). **Fails on the other zero-weight states** (e.g. states with tau_F = +-i/4 have W -> 0 and K_1 = 0 but the shortcut gives B_1 = -+ i), so it is not claimed beyond the singlets.
* General formula test [hp 90-100 digits]: with tau' the eigenvalue of t_1'(theta_F1) on the eps-chain, |W - B_1^2|/(1 + |B_1^2|) over **all 28 zero-weight states**: max 3.3e-3 at eps = 1/800, 1.6e-3 at eps = 1/1600 (halves with eps), i.e. W -> B_1^2 with B_1 = 1 + 4 tau + 2 tau' + 16 tau^2 on every state.

Final L = 2 formula with every piece at the untwisted point:
|C_{M=1}| = K_1 B_1^2 / sqrt(B_0), K_1 = 2 - A - B - C (one-pair Q-data), B_1 = 1 + 4 tau_1(theta_F) + 2 tau_1'(theta_F) + 16 tau_1(theta_F)^2, B_0 = lim B the known SoV norm (2.5e6, 7776 on the singlets).
(Check: 10 * 25 / sqrt(2.5e6) = 1/sqrt40 ... with B_0 = 2.5e6: 250/1581.14 = 0.158114; 2 * 9 / sqrt(7776) = 0.204124.)

## Not done / caveats
* B (the SoV norm) not limited or identified (scope). tau_1'(theta_F) as an expression in Q's (Q_j'/Q_j at theta_F) not worked out; the parity shortcut removes it only on the singlets.
* General-L proof of B_z; L = 3 norm; middle bridges; determinant forms: open. K_kappa at z2 != z3 holds only in the limit (the control residual is 0.069 at generic z).
* The finite-twist C_Q equals C_direct exactly (eigenvector norm), and approaches the data only as eps -> 0 (regularised structure constant, as in the sources); with the Q-norm B they differ by sqrt(H), H = 1 + O(eps).
* The Gram matrix built here (diagonal frame twist, companion-frame words) has 5694 non-zeros at eps = 1/50, not the 3323 quoted in `agent-2026-10-13-sov-norm.md`; B = sum tauT M tau still equals the direct pairing to 1e-59 on all 28 states, so the basis is a valid (different) regularised basis or the count there refers to another normalisation of the words. Not investigated (B is out of scope).
* Kernel memory: the sparse L = 2 monodromy builders retain memory after use (about 0.1 GB per eps-build set, not released by Share/ClearSystemCache/clearing my symbols; my own symbols are < 10 MB); the shared kernel went from 0.07 to about 4.5 GB over the session. No kernel was restarted, aborted or created.
