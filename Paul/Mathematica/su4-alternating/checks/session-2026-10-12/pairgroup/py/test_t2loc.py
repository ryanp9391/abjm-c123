# prediction: at theta_A - theta_F = 2h (pair 1), t2(theta_F + h) localises on pair 1 for any global rotation
import numpy as np
from pgbase import *
rng = np.random.default_rng(7)
h = 0.5
z = np.array([1.5, 1.25, 2/3, 4/7])
def bpspair(v, u):  # project u so that v.u = 0
    u = u - (v@u)/(v@v)*v; return v, u
def resid_local(out, rest, npair=1):
    """distance of out from (anything on first 2*npair sites) (x) rest"""
    M = out.reshape(4**(2*npair), -1)
    r = rest.reshape(-1)
    proj = np.outer(M @ r / (r@r), r)
    return np.linalg.norm(M - proj)/np.linalg.norm(M)
for label, th in (("ABJM pair1", [-0.6+2*h, -0.6, 8/11, -11/8]), ("generic", [2/7, -0.6, 8/11, -11/8])):
    for G in (np.diag(z), rng.normal(size=(4,4))):
        ch = Chain(['A','F','A','F'], th, h, G.astype(object))
        v, u = bpspair(rng.normal(size=4), rng.normal(size=4))
        psi = prod_state([v, u, v, u], dtype=float)
        rest = prod_state([v, u], dtype=float)
        for k, sh, base in ((1,0,1),(2,h,1),(2,0,1),(2,2*h,1),(3,-h,0),(3,0,1),(3,h,1),(3,2*h,1),(3,3*h,1),(2,-h,0),(1,-h,0),(1,h,1)):
            uu = th[base] + sh
            out = ch.apply_t(k, uu, psi)
            print(label, 'G diag' if np.allclose(G, np.diag(np.diag(G))) else 'G rand', 't%d(theta%d%+.1f)'%(k, base+1, sh), 'nonlocal resid %.2e'%resid_local(out, rest))
        print()
