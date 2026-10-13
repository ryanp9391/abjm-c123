# Task 4d: Q-form of the pair-1 eigenvalues at L = 3 pairs (ABJM point, distinct theta_F, generic twist), states of the
# weight sectors (2,0,0,-2) [24 states] and (3,-1,-1,-1)... using sector matrices; Q's from the Baxter equation.
import numpy as np, time
import eng, lp, qf, qdict, untw
h = 0.5; thF = [-0.6, -1.375, -7 / 9]; z = [1.5, 1.25, 2 / 3, 4 / 7]
ch = lp.chain_abjm(3, thF, h, z, 'c')
us = np.linspace(-2.3, 2.1, 34) + 0.0123
for w in [(2, 0, 0, -2), (2, 1, -1, -2)]:
    t0 = time.time()
    sec = untw.SECT[w]
    Mg = untw.sector_matrix(ch, untw.GEN, sec)
    ev, V = np.linalg.eig(Mg)
    def eig(k, u, reg=False):
        M = untw.sector_matrix(ch, [(1.0, k, u, reg)], sec); MV = M @ V
        return np.einsum('ij,ij->j', V.conj(), MV) / np.einsum('ij,ij->j', V.conj(), V)
    taus = {0: [np.ones(len(sec), complex)] * len(us)}
    for k in (1, 2, 3, 4): taus[k] = [eig(k, u) for u in us]
    uF = thF[0]; u0 = uF + h
    tF, x, t3h = eig(1, uF), eig(1, u0), eig(3, u0, True); t2 = eig(2, u0)
    worst = np.zeros(4); fails = 0
    for s in range(len(sec)):
        ta = [[taus[a][n][s] for n in range(len(us))] for a in range(5)]
        qs = [qf.solve_q(ch, ta, us, z[i], 10)[0] for i in range(4)]
        if any(q is None for q in qs): fails += 1; continue
        for I in [(0, 1, 2), (1, 2, 3)]:
            fF, fx, f2, f3, res = qdict.formulas(ch, qs, z, 1, I)
            rel = lambda a, b: abs(a - b) / max(abs(b), 1e-12)
            worst = np.maximum(worst, [max(rel(t, tF[s]) for t in fF), rel(fx, x[s]), rel(f2, t2[s]), rel(f3, t3h[s])])
    print('L=3 sector %s (%d states, %d Baxter failures): max rel. dev. tau_F %.1e, x %.1e, tau2(u0) %.1e, tau3^(u0) %.1e  (%.0f s)'
          % (w, len(sec), fails, *worst, time.time() - t0))
