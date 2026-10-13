# for each pair geometry d = theta_A - theta_F (in units of h), which single t_k(theta_F + m h) localise on pair 1
# (4 sites, random twist, random BPS global rotation); binary-exact floats so reg works
import numpy as np, itertools
from pgbase import *
rng = np.random.default_rng(5)
h = 0.5; thF = -0.625; th3, th4 = 0.8125, -1.3125
G = rng.normal(size=(4,4))
def bps():
    v = rng.normal(size=4); u = rng.normal(size=4); u = u - (v@u)/(v@v)*v; return v, u
v, u = bps()
ref = prod_state([v,u,v,u], dtype=float); rest = np.kron(v,u)
def nonloc(out):
    M = out.reshape(16,16); n = np.linalg.norm(M)
    if n < 1e-12: return -1
    return np.linalg.norm(M - np.outer(M@rest/(rest@rest), rest))/n
for d in range(-5, 6):
    thA = thF + d*h
    ch = Chain(['A','F','A','F'], [thA, thF, th3, th4], h, G.astype(object))
    found = []
    for k in (1,2,3):
        for m in range(-6, 8):
            out = ch.apply_t(k, thF + m*h, ref, reg=True)
            r = nonloc(out)
            if r == -1: found.append('t%d(F%+d)=0' % (k, m))
            elif r < 1e-9: found.append('t%d(F%+d)' % (k, m))
    print('d=%+d:' % d, ', '.join(found))
