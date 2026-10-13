import time, itertools
import numpy as np
from fractions import Fraction as Fr
from sov_l4 import *

# float chain
hf = 0.5; thf = tuple(float(x) for x in TH0); zf = tuple(float(x) for x in Z0)
chf = Chain(h=hf, th=thf, z=zf, exact=False)
t0 = time.time()
A = chf.t(1, 0.3); B = chf.t(3, -0.7); C = chf.t(2, 0.45)
print('build float t1,t2,t3', time.time() - t0)
print('[t1,t3]', np.abs(A @ B - B @ A).max(), '[t1,t2]', np.abs(A @ C - C @ A).max(), '[t2,t3]', np.abs(B @ C - C @ B).max())

# quantum minors check (float)
def Tmats(u):
    P = monodromy_dicts(u, thf, hf)
    T = [[to_dense(P[i][j], complex) for j in range(4)] for i in range(4)]
    G = np.array(Gc(zf), dtype=complex)
    TG = [[sum(T[i][k] * G[k, j] for k in range(4)) for j in range(4)] for i in range(4)]
    return TG

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s

def ta_qmin(a, u):
    Ts = [Tmats(u - k * hf) for k in range(a)]
    tot = np.zeros((256, 256), complex)
    for I in itertools.combinations(range(4), a):
        for p in itertools.permutations(range(a)):
            M = np.identity(256, complex)
            for k in range(a):
                M = M @ Ts[k][I[p[k]]][I[k]]
            tot += perm_sign(p) * M
    return tot

for a, u in [(2, 0.45), (3, -0.7)]:
    q = ta_qmin(a, u)
    f = chf.t(a, u)
    print('fused vs qmin a=%d' % a, np.abs(q - f).max(), np.abs(f).max())
q4 = ta_qmin(4, 0.31)
print('t4 scalar?', np.abs(q4 - q4[0, 0] * np.identity(256)).max(), q4[0, 0])
