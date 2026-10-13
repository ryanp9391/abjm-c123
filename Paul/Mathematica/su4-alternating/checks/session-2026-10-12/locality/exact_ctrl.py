import sys, time
from fractions import Fraction as Fr
import numpy as np
from lib import *
L = int(sys.argv[1]); thp = int(sys.argv[2]); a0 = Fr(1, 3)
th = {1: TH6, 2: TH6b, 3: TH6c}[thp]
ch = Chain(L, th, H, np.diag(np.array(Z4, dtype=object)) if False else None, Fraction) if False else build_chain(L, which_th=1, dt=Fraction)
ch.th = th[:L]
h = H; z = Z4
def kron_obj(vs):
    out = np.array([Fr(1)], dtype=object)
    for v in vs:
        out = np.array([x * y for x in out for y in v], dtype=object)
    return out
v = [Fr(1), Fr(0), Fr(0), a0]; u = [-a0, Fr(0), Fr(0), Fr(1)]
def W_exact(L, first=None):
    vs = []
    for p in range(L // 2):
        if p == 0 and first is not None: vs += first
        else: vs += [v, u]
    return kron_obj(vs)
psi0 = W_exact(L)
def site_scale(psi, site, d):
    """psi with site-vector components multiplied by d (list of 4) (acting as psi.rho(diag d) on that site)"""
    n = L
    arr = psi.reshape(4 ** site, 4, 4 ** (n - site - 1)).copy()
    for i in range(4): arr[:, i, :] = arr[:, i, :] * d[i]
    return arr.reshape(-1)
def vec_eq_scalar(x, y):
    """is x = c*y ? return c or None"""
    nz = [i for i in range(len(y)) if y[i] != 0]
    c = x[nz[0]] / y[nz[0]]
    for i in range(len(x)):
        if x[i] != c * y[i]: return None
    return c
t0 = time.time()
# c1 formula
c1inv = h
for s in range(1, L):
    c1inv *= (th[0] - th[s]) if ch.types[s] == 'F' else (th[0] - th[s] + h)
c1 = 1 / c1inv
lhs = ch.apply(1, th[0], site_scale(psi0, 0, [1 / zz for zz in z]) * c1)
print('calibration  psi0 = c1 psi0 rho1(Lam^-1) t1(th1):', vec_eq_scalar(lhs, psi0) == 1, 'time', time.time() - t0, flush=True)
# t1(th1)^-1 t1(th2) -> c12 rho1(L^-1) rho2(L)
b = site_scale(psi0, 0, [1 / zz for zz in z]) * c1      # psi0 t1(th1)^-1
bb = ch.apply(1, th[1], b)
tgt = site_scale(site_scale(psi0, 0, [1 / zz for zz in z]), 1, list(z))
c12 = vec_eq_scalar(bb, tgt)
c12f = Fr(1)
for s in range(2, L):
    c12f *= (th[1] - th[s]) / (th[0] - th[s] + h) if False else ((th[1] - th[s]) / (th[0] - th[s] + h) if ch.types[s] == 'A' else (th[1] - th[s] - h) / (th[0] - th[s]))
print('c12 exact', c12, ' conj formula', c12f, ' equal:', c12 == c12f, flush=True)
# t3(th1-h) on psi0
y = ch.apply(3, th[0] - h, psi0)
# is y = c * psi0 rho1(diag d)?  find d by ratios of components e1, e4 sites
arr_y = y.reshape(4, -1); arr_0 = psi0.reshape(4, -1)
ratios = []
for i in range(4):
    nz = [k for k in range(arr_0.shape[1]) if arr_0[i, k] != 0]
    if not nz: ratios.append(None); continue
    r = arr_y[i, nz[0]] / arr_0[i, nz[0]]
    ok = all(arr_y[i, k] == r * arr_0[i, k] for k in range(arr_0.shape[1]))
    ratios.append((r, ok))
print('t3(th1-h) psi0 = psi0 rho1(diag(mu_i)): mu =', ratios)
if ratios[0] and ratios[3]:
    mu1, mu4 = ratios[0][0], ratios[3][0]
    print(' mu4/mu1 =', mu4 / mu1, ' z1/z4 =', z[0] / z[3], ' z4/z1 =', z[3] / z[0])
print('time', time.time() - t0)
