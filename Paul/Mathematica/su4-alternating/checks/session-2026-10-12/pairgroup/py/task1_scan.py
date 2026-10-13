# Task 1: on the good companion reference <0| = e1 e4 e1 e4 (4 sites), which t_k(theta_s + m h) localise on a single site?
import numpy as np, itertools, sys
from pgbase import *
h = Fr(1,2); th = [Fr(2,7), Fr(-3,5), Fr(8,11), Fr(-11,8)]
if len(sys.argv) > 1 and sys.argv[1] == 'abjm':
    th = [Fr(-3,5)+1, Fr(-3,5), Fr(8,11)+1, Fr(8,11)]
z = [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
Gc = companion(z)
ch = Chain(['A','F','A','F'], th, h, Gc)
e = lambda i: [Fr(int(j == i)) for j in range(4)]
vac = [e(0), e(3), e(0), e(3)]
ref = np.vectorize(Fr, otypes=[object])(prod_state(vac))
def site_loc(out):
    """if out = (vector on site s) (x) vac elsewhere for a single s, return (s, vector)"""
    res = []
    for s in range(4):
        others = [vac[t] for t in range(4) if t != s]
        o = np.moveaxis(out, s, 0).reshape(4, -1)
        r = prod_state(others).reshape(-1); r = np.vectorize(Fr, otypes=[object])(r)
        vec = o.dot(r)  # r is unit (basis vectors)
        rec = np.multiply.outer(vec, r)
        if all(x == 0 for x in (o - rec).reshape(-1)) and any(x != 0 for x in vec):
            res.append((s, vec))
    return res
for k in (1,2,3):
    for s in range(4):
        for m in range(-3, 5):
            u0 = th[s] + m*h
            out = ch.apply_t(k, u0, ref, reg=True)
            if all(x == 0 for x in out.reshape(-1)): continue
            loc = site_loc(out)
            if all(x == 0 for x in (out - ref * (out.reshape(-1)[np.argmax([x != 0 for x in ref.reshape(-1)])])).reshape(-1)):
                continue  # eigen (scalar)
            if loc:
                for (ss, vec) in loc:
                    print('t%d(theta%d%+dh): site %d -> %s' % (k, s+1, m, ss+1, [str(x) for x in vec]))
