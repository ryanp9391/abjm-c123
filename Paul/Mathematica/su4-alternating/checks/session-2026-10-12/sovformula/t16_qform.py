"""Task 4: Q-form of the localised-word eigenvalues at the ABJM point (pair by pair), L = 2, generic diagonal twist.
Gauge (sov_l4): Q_i(u) = z_i^{u/h} Phi(u) q_i(u), Phi(u) = V4(u) Phi(u-h), V4(u) = prod_A (u-thA) prod_F (u-thF-h);
Baxter: sum_a (-1)^a z_i^{-a} tau_a(u) q_i(u-ah) Pi_a(u) = 0, Pi_a(u) = prod_{b=a}^{3} V4(u-bh), Pi_4 = 1.
At u0 = thF1 + h = thA1 - h: tau3(u0) = tau4(u0) = tau4'(u0) = 0, Pi_0(u0) = 0, tau2(u0) = x tauF, and the u-derivative gives
   A_i w + B_i x + C_i y + K_i = 0   (i = 1..4, rank 3), with
   A_i = z_i^-2 q_i(u0-2h) Pi_2(u0),  B_i = -z_i^-1 [q_i(u-h) Pi_1(u)]'(u0) + z_i^-2 tauF [q_i(u-2h) Pi_2(u)]'(u0),
   C_i = -z_i^-3 q_i(u0-3h) Pi_3(u0), K_i = q_i(u0) V4'(u0) Pi_1(u0),
   x = tau1(u0), w = tau2'(u0) - tau1'(u0) tauF, y = tau3'(u0) (= eigenvalue of Z = t3hat(u0)), tauF = tau1(thF1) = z_i q_i(thF1) V4(thF1)/q_i(thF1-h)."""
import numpy as np, itertools, mpmath as mp, sys
from fractions import Fraction as Fr
import sov_l4 as S
mp.mp.dps = 50
h = Fr(1, 2); thF = (Fr(-3, 5), Fr(-11, 8)); th = (thF[0] + 2 * h, thF[0], thF[1] + 2 * h, thF[1])
z = S.Z0
chd = S.Chain(h=h, th=th, G=S.Gdiag(z))
tm = S.tomp
H = tm(h); TH = [tm(t) for t in th]; Z = [tm(x) for x in z]
UP = [Fr(k, 3) + Fr(1, 17) for k in range(-2, 3)]
def V4(u): return mp.fprod([(u - TH[s]) if S.TYPES[s] == 'A' else (u - TH[s] - H) for s in range(4)])
def g2(u): return mp.fprod([(u - TH[s] - H) if S.TYPES[s] == 'A' else (u - TH[s]) for s in range(4)])
def g3(u): return mp.fprod([(u - TH[s] - H) * (u - TH[s] - 2 * H) if S.TYPES[s] == 'A' else (u - TH[s]) * (u - TH[s] - H) for s in range(4)])
def qd(s, u):
    t = TH[s]
    return (u - t + H) * (u - t - H) * (u - t - 2 * H) * (u - t - 3 * H) if S.TYPES[s] == 'A' else (u - t) * (u - t - H) * (u - t - 2 * H) * (u - t - 4 * H)
def peval(c, u):
    r = mp.mpf(0)
    for x in reversed(c): r = r * u + x
    return r
blocks = S.weight_blocks()
mats = {(a, u): chd.t(a, u, reduced=True) for a in (1, 2, 3) for u in UP}
u0 = th[1] + h
T1u0 = chd.t(1, u0); T1F = chd.t(1, th[1]); T3hat = chd.t(3, u0, reduced=True)   # t3/g3 ; t3hat = (t3/g3) * g3 without the vanishing factor
g3hat = mp.fprod([(tm(u0) - TH[s] - H) * (tm(u0) - TH[s] - 2 * H) if S.TYPES[s] == 'A' else ((tm(u0) - TH[s]) * (tm(u0) - TH[s] - H) if s != 1 else (tm(u0) - TH[s])) for s in range(4)])
Mc = chd.t(1, Fr(3, 11)) + Fr(7, 10) * chd.t(3, Fr(-5, 7), reduced=True) - Fr(3, 13) * chd.t(2, Fr(2, 9), reduced=True)
Vm = mp.matrix([[tm(u) ** j for j in range(5)] for u in UP])
def solve_q(tau, i, Dmax=8):
    zi = Z[i]
    samp = [mp.mpf(k) / 7 + mp.mpf(1) / 29 for k in range(-20, 21)]
    for Dg in range(Dmax + 1):
        rows = []
        for u in samp[:22 + Dg]:
            row, sc = [mp.mpf(0)] * (Dg + 1), mp.mpf(0)
            for a in range(5):
                pre = (-1) ** a * zi ** (-a) * tau(a, u) * mp.fprod([V4(u - b * H) for b in range(a, 4)])
                x = u - a * H
                for k in range(Dg + 1): row[k] += pre * x ** k
                sc = max(sc, abs(pre) * max(abs(x) ** k for k in range(Dg + 1)))
            rows.append([r / sc for r in row])
        U, Sv, Vt = mp.svd_r(mp.matrix(rows))
        if Sv[Sv.rows - 1] < mp.mpf(10) ** -30:
            v = [Vt[Vt.rows - 1, c] for c in range(Vt.cols)]
            return [x / v[-1] for x in v]
    raise RuntimeError
def dpeval(c, u):
    return peval([k * c[k] for k in range(1, len(c))], u) if len(c) > 1 else mp.mpf(0)
U0 = tm(u0); TF1 = TH[1]
def dnum(f, u, eps=mp.mpf(10) ** -20): return (f(u + eps) - f(u - eps)) / (2 * eps)
maxerr = {'tauF': 0, 'x': 0, 'w': 0, 'y': 0}
nstates = 0
for wgt, idx in sorted(blocks.items()):
    if len(idx) > 6 or nstates >= 40: continue
    B = mp.matrix([[tm(Mc[r, c]) for c in idx] for r in idx])
    E, ER = mp.eig(B)
    for k in range(len(idx)):
        v = ER[:, k]; p = max(range(len(idx)), key=lambda r: abs(v[r])); v = v / v[p]
        def ev(M):
            tv = mp.matrix([[tm(M[r, c]) for c in idx] for r in idx]) * v
            return tv[p] / v[p]
        tb = {}
        for a in (1, 2, 3):
            vals = [ev(mats[(a, u)]) for u in UP]
            c = mp.lu_solve(Vm, mp.matrix(vals)); tb[a] = [c[j] for j in range(5)]
        def tau(a, u):
            if a == 0: return mp.mpf(1)
            if a == 1: return peval(tb[1], u)
            if a == 2: return g2(u) * peval(tb[2], u)
            if a == 3: return g3(u) * peval(tb[3], u)
            return mp.fprod(Z) * mp.fprod([qd(s, u) for s in range(4)])
        # direct eigenvalues
        x_dir = ev(T1u0); tauF_dir = ev(T1F); y_dir = ev(T3hat) * g3hat
        tau1p = dnum(lambda u: tau(1, u), U0); tau2p = dnum(lambda u: tau(2, u), U0)
        w_dir = tau2p - tau1p * tauF_dir
        q = [solve_q(tau, i) for i in range(4)]
        Pi = lambda a, u: mp.fprod([V4(u - b * H) for b in range(a, 4)])
        rows = []
        for i in range(4):
            qi = lambda u, i=i: peval(q[i], u)
            A = Z[i] ** -2 * qi(U0 - 2 * H) * Pi(2, U0)
            Bc = -Z[i] ** -1 * dnum(lambda u: qi(u - H) * Pi(1, u), U0) + Z[i] ** -2 * tauF_dir * dnum(lambda u: qi(u - 2 * H) * Pi(2, u), U0)
            C = -Z[i] ** -3 * qi(U0 - 3 * H) * Pi(3, U0)
            K = qi(U0) * dnum(V4, U0) * Pi(1, U0)
            rows.append((A, Bc, C, K))
        sol = mp.lu_solve(mp.matrix([r[:3] for r in rows[:3]]), mp.matrix([-r[3] for r in rows[:3]]))
        w_Q, x_Q, y_Q = sol[0], sol[1], sol[2]
        res4 = abs(rows[3][0] * w_Q + rows[3][1] * x_Q + rows[3][2] * y_Q + rows[3][3]) / max(abs(rows[3][3]), mp.mpf(1) ** 0)
        tauF_Q = [Z[i] * peval(q[i], TF1) * V4(TF1) / peval(q[i], TF1 - H) for i in range(4)]
        rel = lambda a, b: abs(a - b) / max(abs(b), mp.mpf(10) ** -30)
        maxerr['tauF'] = max(maxerr['tauF'], max(rel(t, tauF_dir) for t in tauF_Q))
        maxerr['x'] = max(maxerr['x'], rel(x_Q, x_dir)); maxerr['w'] = max(maxerr['w'], rel(w_Q, w_dir)); maxerr['y'] = max(maxerr['y'], rel(y_Q, y_dir))
        nstates += 1
        print(wgt, 'x', mp.nstr(x_dir, 8), 'xQ', mp.nstr(x_Q, 8), '| y', mp.nstr(y_dir, 8), 'yQ', mp.nstr(y_Q, 8), '| w', mp.nstr(w_dir, 8), 'wQ', mp.nstr(w_Q, 8), '| 4th eq res', mp.nstr(res4, 3), flush=True)
print('states:', nstates, ' max rel errors:', {k: mp.nstr(v, 3) for k, v in maxerr.items()})
