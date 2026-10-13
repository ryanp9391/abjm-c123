"""Twisted Q-functions from the 4th-order Baxter equation (mpmath).
Baxter operator  O = sum_{a=0}^4 (-1)^a tau_a(u) D^a,  D f(u) = f(u-h)  (tau_a = full fused eigenvalues, tau_4 = det(Lam) qdet).
On the vacuum O = (1 - z1 V1 D)(1 - z2 V2 D)(1 - z3 V3 D)(1 - z4 V4 D), V_i = vacuum eigenvalues of T_ii on the ket
vacuum e1 e4 e1 e4 (diag frame).  Single-index Q's: O Q_i = 0, Q_i(u) = z_i^{u/h} Phi4(u) q_i(u), Phi4(u) = V4(u) Phi4(u-h),
q_i polynomial; the polynomial equation is
    sum_a (-1)^a z_i^{-a} tau_a(u) q_i(u - a h) prod_{b=a}^{3} V4(u - b h) = 0.
"""
import pickle, itertools, math
import mpmath as mp
from fractions import Fraction as Fr
from sov_l4 import H0, TH0, Z0, TYPES

mp.mp.dps = 50


def tomp(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x)


h = tomp(H0)
th = [tomp(t) for t in TH0]
z = [tomp(t) for t in Z0]
thA = [th[0], th[2]]
thF = [th[1], th[3]]


def peval(c, u):
    r = mp.mpf(0)
    for x in reversed(c):
        r = r * u + x
    return r


def V(i, u):
    """vacuum eigenvalues of T_ii (i = 1..4) on the ket vacuum, generic theta per site"""
    pa = [(u - t) for t in thA]
    pf = [(u - t) for t in thF]
    if i == 1:
        return mp.fprod([x + h for x in pa]) * mp.fprod(pf)
    if i in (2, 3):
        return mp.fprod(pa) * mp.fprod(pf)
    return mp.fprod(pa) * mp.fprod([x - h for x in pf])


def g2(u):
    return mp.fprod([(u - t - h) for t in thA]) * mp.fprod([(u - t) for t in thF])


def g3(u):
    return mp.fprod([(u - t - h) * (u - t - 2 * h) for t in thA]) * mp.fprod([(u - t) * (u - t - h) for t in thF])


def qdetA(u, t):
    return (u - t + h) * (u - t - h) * (u - t - 2 * h) * (u - t - 3 * h)


def qdetF(u, t):  # fitted/verified in check_qdet()
    return (u - t + h) * (u - t - h) * (u - t - 2 * h) * (u - t - 3 * h) * 0 + (u - t) * (u - t - h) * (u - t - 2 * h) * (u - t - 4 * h)


def tau4(u):
    return mp.fprod(z) * mp.fprod([qdetA(u, t) for t in thA]) * mp.fprod([qdetF(u, t) for t in thF])


class State:
    def __init__(self, d):
        self.d = d
        self.weight = d['weight']
        self.tb = d['taubar']

    def tau(self, a, u):
        if a == 0:
            return mp.mpf(1)
        if a == 1:
            return peval(self.tb[1], u)
        if a == 2:
            return g2(u) * peval(self.tb[2], u)
        if a == 3:
            return g3(u) * peval(self.tb[3], u)
        if a == 4:
            return tau4(u)


def baxter_row(st, zi, u, D, with_scale=False):
    """coefficients (in q's monomials c_0..c_D) of sum_a (-1)^a zi^{-a} tau_a(u) q(u-ah) prod_{b=a}^3 V4(u-bh);
    with_scale: also return the size of the individual terms (for a scale-aware rank test)"""
    row = [mp.mpf(0)] * (D + 1)
    scale = mp.mpf(0)
    for a in range(5):
        pre = (-1) ** a * zi ** (-a) * st.tau(a, u) * mp.fprod([V(4, u - b * h) for b in range(a, 4)])
        x = u - a * h
        for k in range(D + 1):
            row[k] += pre * x ** k
        scale = max(scale, abs(pre) * max(abs(x) ** k for k in range(D + 1)))
    return (row, scale) if with_scale else row


SAMPLE = [mp.mpf(k) / 7 + mp.mpf(1) / 29 for k in range(-20, 21)]


def solve_q(st, zi, Dmax=10, tol=mp.mpf(10) ** -35):
    """minimal-degree polynomial solution q (monic, coefficient list low->high) of the twisted Baxter equation.
    Rows are scaled by the size of the individual terms, so a small singular value means genuine cancellation."""
    for D in range(0, Dmax + 1):
        npts = 16 + D + 6
        rows = []
        for m in range(npts):
            r, s = baxter_row(st, zi, SAMPLE[m], D, with_scale=True)
            rows.append([x / s for x in r])
        A = mp.matrix(rows)
        U, S, Vt = mp.svd_c(A)
        smin = S[S.rows - 1]
        if smin < tol:
            v = [mp.conj(Vt[Vt.rows - 1, c]) for c in range(Vt.cols)]
            lead = v[-1]
            v = [x / lead for x in v]
            return v, smin
    return None, None


def check_qdet(states):
    """tau4 vs the 4x4 Baxter consistency: on the vacuum, tau_a = e_a(y_i(u - shifts))"""
    vac = [s for s in states if s.weight == (2, 0, 0, -2)][0]
    u = mp.mpf('0.3141')
    y = lambda i, x: z[i - 1] * V(i, x)
    e = {1: sum(y(i, u) for i in range(1, 5)),
         2: sum(y(i, u) * y(j, u - h) for i, j in itertools.combinations(range(1, 5), 2)),
         3: sum(y(i, u) * y(j, u - h) * y(k, u - 2 * h) for i, j, k in itertools.combinations(range(1, 5), 3)),
         4: y(1, u) * y(2, u - h) * y(3, u - 2 * h) * y(4, u - 3 * h)}
    return {a: mp.nstr(abs(vac.tau(a, u) - e[a]) / abs(e[a]), 5) for a in (1, 2, 3, 4)}


if __name__ == '__main__':
    data = pickle.load(open('eigen.pkl', 'rb'))
    states = [State(d) for d in data['states']]
    print('vacuum check tau_a = e_a(y):', check_qdet(states))
