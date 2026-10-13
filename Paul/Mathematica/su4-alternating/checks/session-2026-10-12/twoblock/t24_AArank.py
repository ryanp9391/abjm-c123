"""Rank of the A-site coupling matrix C(r2,r4)_{a1,a3} = sum_{k2,k4} c_{a1 k2 a3 k4} r2^k2 r4^k4 at random r (generic rank over Q(r2,r4)),
and of the F-site coupling F(tau_A1, tau_A3)_{k2,k4}."""
import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; R = X.R; rng = np.random.default_rng(2)
z = z0
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
vz4 = X.vm(e(4), X.VT); l1 = X.vm(e(1), X.Vi); lam1 = X.r(lam[0])
e1, e4 = np.array(e(1)), np.array(e(4))
def pairb(b):
    b = X.r(b); return [(e1 + b * lam1 % p * vz4) % p, (e4 - b * l1) % p]
def rnd(): return int(rng.integers(1, p))
def bps():
    v = rng.integers(1, p, 4); u = rng.integers(1, p, 4); s = int((v[:3] * u[:3] % p).sum() % p); u[3] = (-s * pow(int(v[3]), -1, p)) % p; return [v, u]
def generic_ranks(c):
    c4 = c.reshape(4, 4, 4, 4)
    out = []
    for _ in range(3):
        r2, r4 = rnd(), rnd()
        pw2 = np.array([pow(r2, k, p) for k in range(4)]); pw4 = np.array([pow(r4, k, p) for k in range(4)])
        Cm = np.einsum('akbl,k,l->ab', c4 % p, pw2, pw4) % p if False else None
        Cm = np.zeros((4, 4), dtype=np.int64)
        for a1 in range(4):
            for a3 in range(4):
                s = 0
                for k2 in range(4):
                    for k4 in range(4):
                        s = (s + int(c4[a1, k2, a3, k4]) * int(pw2[k2]) % p * int(pw4[k4])) % p
                Cm[a1, a3] = s
        t1, t3 = [rnd() for _ in range(4)], [rnd() for _ in range(4)]
        Fm = np.zeros((4, 4), dtype=np.int64)
        for k2 in range(4):
            for k4 in range(4):
                s = 0
                for a1 in range(4):
                    for a3 in range(4):
                        s = (s + int(c4[a1, k2, a3, k4]) * t1[a1] % p * t3[a3]) % p
                Fm[k2, k4] = s
        out.append((modrank(Cm, p), modrank(Fm, p)))
    return out
for name, vs in [('physical b (block1), vac', pairb(Fr(3, 7)) + [e1, e4]),
                 ('two-block b, b\'', pairb(Fr(3, 7)) + pairb(Fr(-2, 5))),
                 ('global b', pairb(Fr(3, 7)) + pairb(Fr(3, 7))),
                 ('generic two-block BPS', bps() + bps())]:
    c = X.coeffs(X.bra(vs))
    print(f'{name:30s} generic rank of A-coupling C(r)_(a1,a3), F-coupling F(tauA)_(k2,k4):', generic_ranks(c))
