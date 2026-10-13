import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; rng = np.random.default_rng(11)
def rv(n): return np.array([X.rnd(rng) for _ in range(n)], dtype=np.int64)
def bps():
    v = rv(4); u = rv(4); s = int((v[:3] * u[:3] % p).sum() % p); u[3] = (-s * pow(int(v[3]), -1, p)) % p; return v, u
e1, e4 = np.array(e(1)), np.array(e(4))
v, u = bps(); v2, u2 = bps()
def show(name, vs):
    c = X.coeffs(X.bra(vs)); print(f'{name:50s}', X.ranks(c))
show('block {1,2} generic, rest vac', [v, u, e1, e4])
show('block {4,1} generic (v on 1, u on 4), rest vac', [v, e4, e1, u])
show('block {2,3} generic (u on 2, v on 3), rest vac', [e1, u, v, e4])
show('block {3,4} generic, rest vac', [e1, e4, v, u])
show('blocks {4,1} and {2,3} generic', [v, u2, v2, u])
show('blocks {1,2} and {3,4} generic', [v, u, v2, u2])
show('site 1 generic v, rest vac', [v, e4, e1, e4])
show('site 4 generic u, rest vac', [e1, e4, e1, u])
show('site 2 generic u, rest vac', [e1, u, e1, e4])
show('site 3 generic v, rest vac', [e1, e4, v, e4])
