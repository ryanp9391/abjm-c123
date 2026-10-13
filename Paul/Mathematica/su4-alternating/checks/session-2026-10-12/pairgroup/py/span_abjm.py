# span of all pair-1-localising words (length <= 3) at the ABJM point; is the common group rotation in it?
import numpy as np, itertools, sys
from pgbase import *
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
h = 0.5
th = [-0.6 + 2*h, -0.6, 0.7273, 0.7273 - 2*h]
G = rng.normal(size=(4,4)) if (len(sys.argv) < 3 or sys.argv[2] != 'diag') else np.diag([1.5,1.25,2/3,4/7])
ch = Chain(['A','F','A','F'], th, h, G.astype(object))
def bps():
    v = rng.normal(size=4); u = rng.normal(size=4); u = u - (v@u)/(v@v)*v; return v, u
v, u = bps()
if len(sys.argv) > 2 and sys.argv[2] == 'diag':
    a = 1/3; v = np.array([1,0,0,a]); u = np.array([-a,0,0,1.])
psi0 = prod_state([v,u,v,u], dtype=float)
rest = np.kron(v,u)
def nonloc(out):
    M = out.reshape(16, 16); proj = np.outer(M@rest/(rest@rest), rest)
    n = np.linalg.norm(M); return (np.linalg.norm(M-proj)/n if n > 1e-12 else -1.0), M@rest/(rest@rest)
eps = 1e-5
def mk(k, u0, reg):
    if not reg: return lambda p: ch.apply_t(k, u0, p)
    return lambda p: (ch.apply_t(k, u0+eps, p) - ch.apply_t(k, u0-eps, p))/(2*eps)
thF = th[1]
ops = {'t1(F)': mk(1, thF, 0), 't1(F+h)': mk(1, thF+h, 0), 't1(F+2h)': mk(1, thF+2*h, 0), 't1(F-h)': mk(1, thF-h, 0), 't1(F+3h)': mk(1, thF+3*h, 0),
       't2(F+h)': mk(2, thF+h, 0), 't2(F+2h)': mk(2, thF+2*h, 0), 't2(F-h)': mk(2, thF-h, 0), 't2(F+3h)':mk(2,thF+3*h,1), 't2F':mk(2,thF,1),
       't3^(F+h)': mk(3, thF+h, 1), 't3(F+2h)': mk(3, thF+2*h, 0), 't3^(F+3h)': mk(3, thF+3*h, 1), 't3(F-h)': mk(3, thF-h, 0), 't3^F': mk(3,thF,1), 't3(F+4h)':mk(3,thF+4*h,1)}
names = list(ops)
loc = {'1': psi0.reshape(16,16)@rest/(rest@rest)}
cache = {(): psi0}
for L in (1, 2, 3):
    for w in itertools.combinations_with_replacement(range(len(names)), L):
        prev = cache.get(w[:-1])
        if prev is None: continue
        out = ops[names[w[-1]]](prev)
        r, pairstate = nonloc(out)
        cache[w] = out
        if 0 <= r < 1e-6:
            loc['*'.join(names[i] for i in w)] = pairstate
print('number of localising words:', len(loc))
B = np.array(list(loc.values())).T
sv = np.linalg.svd(B, compute_uv=False)
print('rank of span:', np.sum(sv > 1e-7*sv[0]), ' sv:', np.round(sv[:20]/sv[0], 9))
Gi = np.linalg.inv(G)
for Xn, X in (('G',G),('G^T',G.T),('G^-1',Gi),('G^-T',Gi.T)):
    for kap in (0.3, -0.7):
        g = np.eye(4) - kap*X
        tgt = np.kron(v@g, u@np.linalg.inv(g).T)
        c, *_ = np.linalg.lstsq(B, tgt, rcond=None)
        print(Xn, kap, 'resid %.2e' % (np.linalg.norm(B@c-tgt)/np.linalg.norm(tgt)))
for nm in list(loc)[:80]: print(nm)
