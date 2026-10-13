# Calibration: 2603-type chain of 6's (wedge^2 C^4), L=2 sites, quantum minors from full matrices.
# Does 1 - c1 kappa t1(theta1+s1) + c2 kappa^2 t2(theta1+s2) act on <K Z (x) K Z| as the group element wedge2(1-kappa Lam) on site 1,
# for a GENERIC global rotation K?
import numpy as np, itertools
from pgbase import subsets, wedge, DK, perm_sign
rng = np.random.default_rng(4)
h = 0.5; th = [0.3125, -0.6875]
S2 = subsets(2)
D2 = {(a,b): DK[2][(a,b)].astype(float) for a in range(4) for b in range(4)}
def eab(a,b): m = np.zeros((4,4)); m[a,b] = 1; return m
def lax(u, t):  # [i,j] -> 6x6 operator
    L = np.zeros((4,4,6,6))
    for i in range(4):
        L[i,i] += (u-t)*np.eye(6)
        for j in range(4): L[i,j] += h*D2[(i,j)]
    return L
def T(u, G):
    L1 = lax(u, th[0]); L2 = lax(u, th[1])
    TT = np.einsum('ikab,kjcd->ijacbd', L1, L2).reshape(4,4,36,36)
    return np.einsum('ikab,kj->ijab', TT, G)
def tk(u, k, G):
    Ts = [T(u - m*h, G) for m in range(k)]
    tot = np.zeros((36,36))
    for I in itertools.combinations(range(4), k):
        for p in itertools.permutations(range(k)):
            M = np.eye(36)
            for m in range(k): M = M @ Ts[m][I[p[m]], I[m]]
            tot += perm_sign(p)*M
    return tot
z = np.array([1.5, 1.25, 0.625, 0.5625]); Lam = np.diag(z)
K = rng.normal(size=(4,4))
wK = wedge(K, 2, dtype=float).astype(float)
Z = np.zeros(6); Z[0] = 1
KZ = Z @ wK
ref = np.kron(KZ, KZ)
def nonloc(out):
    M = out.reshape(6,6); ps = M@KZ/(KZ@KZ); n = np.linalg.norm(M)
    return np.linalg.norm(M - np.outer(ps, KZ))/n if n > 1e-12 else -1, ps
for k in (1,2,3):
    for m in range(-3, 4):
        r, ps = nonloc(ref @ tk(th[0] + m*h, k, Lam))
        if 0 <= r < 1e-9: print('localises: t%d(theta1%+dh)' % (k, m))
print('--- zeros and regularised ---')
eps = 1e-4
for k in (1,2,3):
    for m in range(-3, 4):
        M0 = tk(th[0] + m*h, k, Lam)
        if np.abs(ref@M0).max() < 1e-10:
            Mr = (tk(th[0]+m*h+eps, k, Lam) - tk(th[0]+m*h-eps, k, Lam))/(2*eps)
            r, ps = nonloc(ref@Mr)
            print('t%d(theta1%+dh) = 0 on ref; derivative nonloc %.1e' % (k, m, r))
print('--- half-integer scan ---')
for k in (1,2,3):
    for m in range(-8, 9):
        u0 = th[0] + m*h/2
        M0 = tk(u0, k, Lam)
        out = ref@M0
        if np.abs(out).max() < 1e-10:
            M0 = (tk(u0+eps, k, Lam) - tk(u0-eps, k, Lam))/(2*eps); out = ref@M0; tag = '(reg)'
        else: tag = ''
        r, ps = nonloc(out)
        if 0 <= r < 1e-6: print('localises: t%d(theta1%+.1fh)%s' % (k, m/2, tag))
print('--- continuous scan of t2 nonlocality ---')
us = np.linspace(th[0]-4*h, th[0]+4*h, 321)
vals = []
for u0 in us:
    out = ref@tk(u0, 2, Lam); n = np.linalg.norm(out)
    r, _ = nonloc(out); vals.append(r)
vals = np.array(vals)
idx = np.argsort(vals)[:8]
for i in idx: print('u-theta1 = %.4f h   nonloc %.3e' % ((us[i]-th[0])/h, vals[i]))
