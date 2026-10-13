import sympy as sp, pickle, mpmath as mp, random
from itertools import combinations
mp.mp.dps = 50
def M(x):
    if isinstance(x, sp.Basic):
        return mp.mpf(x.p) / x.q if x.is_Rational else mp.mpf(sp.N(x, 60))
    if isinstance(x, mp.mpc): return mp.re(x) if abs(mp.im(x)) < mp.mpf(10)**-40 else x
    return x
h = mp.mpf(1) / 2; th1 = mp.mpf(2) / 7; th2 = mp.mpf(-3) / 5
z = [mp.mpf(3) / 2, mp.mpf(5) / 4, mp.mpf(2) / 3, mp.mpf(4) / 7]
th = th1 - th2
V4 = lambda x: (x - th1) * (x - th2 - h)
states = pickle.load(open('qfun.pkl', 'rb'))
for s in states:
    s['ketm'] = [M(x) for x in s['ket']]
    s['taus'] = {a: [M(c) for c in s['tau'][a]] for a in range(1, 5)}
def tau(s, a, x):
    if a == 0: return mp.mpf(1)
    cs = s['taus'][a]; return sum(cs[m] * x**m for m in range(len(cs)))
def qj(s, j, x):
    cs = s['q'][j]; return sum(M(cs[m]) * x**m for m in range(len(cs)))
# Gamma-gauge values, stripped of z_j^{theta/h} and of Gamma(theta1-3h), Gamma(theta2)
def gA(b):  # Gamma(theta1 - b h)/Gamma(theta1 - 3h), b = 0..3 (b can be -1? no)
    return mp.fprod([V4(th1 - l * h) for l in range(b + 1, 4)])
def QA(s, j, b):  # Q_j(theta1 - b h) / (Gamma(theta1-3h) z_j^{theta1/h})
    return gA(b) * z[j]**(-b) * qj(s, j, th1 - b * h)
def gF(m):  # Gamma(theta2 + m h)/Gamma(theta2), m = 0, 1 ; m=-1,-2 via inverse
    if m >= 0: return mp.fprod([V4(th2 + l * h) for l in range(0, m)])
    return 1 / mp.fprod([V4(th2 + l * h) for l in range(m, 0)])
def QF(s, j, m):
    return gF(m) * z[j]**m * qj(s, j, th2 + m * h)
def det(rows): return mp.det(mp.matrix(rows))
