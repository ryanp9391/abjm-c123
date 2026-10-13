import numpy as np, pickle, random, time
from fractions import Fraction
from fastchain import *
d = pickle.load(open('twopair_W.pkl', 'rb')); labs = d['labs']; W = d['W']
t0 = time.time(); Winv = inv_mod(W); print('inverse', time.time() - t0)
print('check W Winv = 1:', np.array_equal(W.dot(Winv) % P, np.eye(256, dtype=object)))
pickle.dump(Winv, open('twopair_Winv.pkl', 'wb'))
def coeffs(v1, u1, v2, u2):
    x = np.kron(np.kron(np.kron(np.array(v1, dtype=object), np.array(u1, dtype=object)), np.array(v2, dtype=object)), np.array(u2, dtype=object)) % P
    c = x.dot(Winv) % P
    C = np.zeros((16, 16), dtype=object)
    for n, (a1, k2, a3, k4) in enumerate(labs):
        C[4 * a1 + k2, 4 * a3 + k4] = c[n]
    return C
random.seed(3)
rnd = lambda: [random.randrange(P) for _ in range(4)]
e = lambda i: [1 if k == i else 0 for k in range(4)]
tests = {
 'generic (v1,u1,v2,u2)': (rnd(), rnd(), rnd(), rnd()),
 'pair 1 = vacuum (e1,e4), pair 2 generic': (e(0), e(3), rnd(), rnd()),
 'pair 1 generic, pair 2 = vacuum': (rnd(), rnd(), e(0), e(3)),
 'pair 1 generic BPS, pair 2 generic BPS': None,
}
v1, u1, v2, u2 = rnd(), rnd(), rnd(), rnd()
# BPS: fix u_4 so that v.u = 0
def bps(v, u):
    u = list(u); u[3] = (-(v[0] * u[0] + v[1] * u[1] + v[2] * u[2]) * pow(v[3], P - 2, P)) % P; return v, u
v1, u1 = bps(v1, u1); v2, u2 = bps(v2, u2)
tests['pair 1 generic BPS, pair 2 generic BPS'] = (v1, u1, v2, u2)
for name, (a, b, c_, d_) in tests.items():
    C = coeffs(a, b, c_, d_)
    print(name, ': rank of 16x16 C[(a1,k2),(a3,k4)] =', rank_mod(C), '; nonzero entries:', int(np.sum(C % P != 0)))
