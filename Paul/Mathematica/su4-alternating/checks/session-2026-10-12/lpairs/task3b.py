# Task 3b (i): word residual for z2 != z3 at fixed twist, L = 3, ABJM point (numeric).
#   K_A: y -> (z2+z3)/2 in the X2 term, y^2 -> z2 z3 in the X3 and X3X1 terms;  also K_B: y = z2;  and least-squares best.
import numpy as np
import eng, lp
h = 0.5; thF = [-0.6, -1.375, -7 / 9]
a = 1 / 3
for kap in (0.3, 1 / (4 / 7)):
    for dp in (1e-1, 1e-2, 1e-3, 1e-4):
        z = [1.5, 1.25 + dp / 2, 1.25 - dp / 2, 4 / 7]
        ch = lp.chain_abjm(3, thF, h, z, 'c')
        ref = lp.Wstate(ch, [a] * 3)
        X = lp.Xs(ch, ref, 1)
        l1, l2, l3 = eng.lam(ch, 1, ch.th[1], 1), eng.lam(ch, 2, ch.th[1] + h, 1), eng.lam(ch, 3, ch.th[1] + h, 1)
        s = 1 - kap * (z[0] + z[3]); ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
        tgt = (1 - kap * z[0]) ** 2 * lp.Wstate(ch, [ap, a, a])
        def K(y1, y2sq):
            return (s * s * ref + kap * s * (X[0] / (h**2 * l1) + X[1] / (2 * h**4 * y1 * l2) + X[2] / (6 * h**5 * y2sq * l3))
                    + kap**2 * X[3] / (24 * h**7 * y2sq * l1 * l3))
        nt = np.linalg.norm(tgt)
        rA = np.linalg.norm(K((z[1] + z[2]) / 2, z[1] * z[2]) - tgt) / nt
        rB = np.linalg.norm(K(z[1], z[1] ** 2) - tgt) / nt
        A = np.stack([v.reshape(-1) for v in [ref] + list(X)], axis=1)
        cfit, *_ = np.linalg.lstsq(A, tgt.reshape(-1), rcond=None)
        rLS = np.linalg.norm(A @ cfit - tgt.reshape(-1)) / nt
        print('kappa=%.4f  z2-z3=%.0e: rel. residual  K_A %.2e   K_B %.2e   least-squares best %.2e' % (kap, dp, rA, rB, rLS))
