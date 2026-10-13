# scan at the ABJM point: which t_k(theta_s + m h) (regularised) and products of two localise on pair 1 (4 sites)
import numpy as np, itertools
from pgbase import *
rng = np.random.default_rng(11)
h = 0.5
th = [-0.6 + 2*h, -0.6, 0.7273, 0.7273 - 2*h]   # both pairs at the ABJM point
G = rng.normal(size=(4,4))
ch = Chain(['A','F','A','F'], th, h, G.astype(object))
def bps():
    v = rng.normal(size=4); u = rng.normal(size=4); u = u - (v@u)/(v@v)*v; return v, u
v, u = bps()
psi0 = prod_state([v,u,v,u], dtype=float)
rest = np.kron(v,u)
def nonloc(out):
    M = out.reshape(16, 16)
    proj = np.outer(M@rest/(rest@rest), rest)
    n = np.linalg.norm(M)
    return np.linalg.norm(M-proj)/n if n > 1e-12 else -1.0
def op(k, u0):
    """returns function psi -> psi.t_k(u0), regularised if identically zero (derivative)"""
    test = ch.apply_t(k, u0, psi0)
    if np.linalg.norm(test) > 1e-9*np.linalg.norm(psi0):
        return lambda p: ch.apply_t(k, u0, p), ''
    eps = 1e-5
    t2 = ch.apply_t(k, u0+eps, psi0)
    if np.linalg.norm(t2) < 1e-9: return None, 'zero2'
    return (lambda p: (ch.apply_t(k, u0+eps, p) - ch.apply_t(k, u0-eps, p))/(2*eps)), 'reg'
ops = {}
for k in (1,2,3):
    for s in (0,1):
        for m in range(-3, 4):
            u0 = th[s] + m*h
            f, tag = op(k, u0)
            if f is None: continue
            ops['t%d(th%d%+d h)%s' % (k, s+1, m, tag)] = f
single = {}
for name, f in ops.items():
    out = f(psi0); r = nonloc(out); single[name] = out
    if r < 1e-6: print('LOCAL single:', name, '%.1e' % r)
print('---- products of two ----')
names = list(ops)
loc2 = []
for a, b in itertools.combinations_with_replacement(names, 2):
    out = ops[b](single[a]); r = nonloc(out)
    if 0 <= r < 1e-6:
        loc2.append((a,b)); print('LOCAL pair:', a, '*', b, '%.1e' % r)
