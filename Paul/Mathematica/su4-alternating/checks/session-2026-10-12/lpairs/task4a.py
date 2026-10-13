# Task 4a: L = 2 pairs, ABJM point, generic diagonal twist: Q-functions of all 256 states and the dictionary for
#   tau_F = tau1(thF), x = tau1(u0), tau2(u0) = x tau_F, tau3^(u0)  (u0 = thF + h, pair 1).
import numpy as np, time, pickle
import eng, lp, qf
h = 0.5; thF = [-0.6, -1.375]
z = [1.5, 1.25, 2 / 3, 4 / 7]
ch = lp.chain_abjm(2, thF, h, z, 'c')
t0 = time.time()
V, ev, gaps = qf.joint_eigen(ch)
print('joint eigenbasis: min gap %.2e (%.0f s)' % (gaps.min(), time.time() - t0))
us = np.linspace(-2.3, 2.1, 30) + 0.0123
taus = {0: [np.ones(256, complex) for _ in us]}
for k in (1, 2, 3, 4):
    taus[k] = []
    for u in us:
        tv, res = qf.eigvals_at(ch, V, k, u)
        taus[k].append(tv)
print('eigenvalues at %d points (%.0f s)' % (len(us), time.time() - t0))
uF = thF[0]; u0 = uF + h
tF, r1 = qf.eigvals_at(ch, V, 1, uF)
xx, r2 = qf.eigvals_at(ch, V, 1, u0)
t2u0, r3 = qf.eigvals_at(ch, V, 2, u0)
t3h, r4 = qf.eigvals_at(ch, V, 3, u0, reg=True)
print('eig residuals: %.1e %.1e %.1e %.1e' % (r1.max(), r2.max(), r3.max(), r4.max()))
print('identity tau2(u0) = x tau_F: max rel dev %.1e' % np.max(np.abs(t2u0 - xx * tF) / np.maximum(np.abs(t2u0), 1e-12)))
# Q's
Q = []
bad = 0
for s in range(256):
    ta = [[taus[a][n][s] for n in range(len(us))] for a in range(5)]
    qs = []
    for i in range(4):
        c, info = qf.solve_q(ch, ta, us, z[i], 8)
        if c is None: bad += 1
        qs.append(c)
    Q.append(qs)
print('Baxter solved: failures %d (%.0f s)' % (bad, time.time() - t0))
pickle.dump(dict(V=V, ev=ev, gaps=gaps, us=us, taus=taus, tF=tF, x=xx, t2u0=t2u0, t3h=t3h, Q=Q, z=z, thF=thF, h=h), open('task4a.pkl', 'wb'))
degs = [tuple(len(q) - 1 for q in qs) for qs in Q if all(q is not None for q in qs)]
from collections import Counter
print('degree patterns:', Counter(degs).most_common(8))
