import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; R = X.R
G = X.Gc; GT = X.GcT; I = R.arr(np.eye(4, dtype=int))
e1, e4 = np.array(e(1)), np.array(e(4))
def lin(k, A): return (I - (X.r(k) * A) % p) % p
def vm(v, A): return X.vm(v, A)
def poly(cs, A):
    out = np.zeros((4, 4), dtype=np.int64); Pk = I.copy()
    for c in cs:
        out = (out + X.r(c) * Pk) % p; Pk = R.mm(Pk, A)
    return out
def flag(dp, dq, sd):
    return [vm(e1, poly([1, Fr(2, 3) + sd, Fr(-5, 4), Fr(3, 7)][:dp + 1], GT)), vm(e4, poly([1, Fr(-3, 5) - sd, Fr(4, 9), Fr(2, 11)][:dq + 1], G))]
v1, u1, v2, u2 = [2, -1, 3, 1], [1, 1, -2, 5], [1, 3, -2, 2], [3, -1, 1, 1]
k1, k2, k3, k4 = Fr(3, 7), Fr(-5, 11), Fr(2, 9), Fr(7, 5)
inv = lambda A: modinv(A, p)
fam = [
 ("generic two-block BPS (v⊗u)⊗(v'⊗u')", [v1, u1, v2, u2]),
 ("generic global BPS (v⊗u)⊗(v⊗u)", [v1, u1, v1, u1]),
 ("block {1,2} generic, block {3,4} = SoV vacuum", [v1, u1, e1, e4]),
 ("block {3,4} generic, block {1,2} = SoV vacuum", [e1, e4, v1, u1]),
 ("block {4,1} generic (A₁ = v, F₄ = u), rest vacuum", [v1, e4, e1, u1]),
 ("block {2,3} generic (F₂ = u, A₃ = v), rest vacuum", [e1, u1, v1, e4]),
 ("site 1: e₁(1−κGᵀ)", [vm(e1, lin(k1, GT)), e4, e1, e4]),
 ("site 2: e₄(1−κG)", [e1, vm(e4, lin(k2, G)), e1, e4]),
 ("site 3: e₁(1−κGᵀ)", [e1, e4, vm(e1, lin(k3, GT)), e4]),
 ("site 4: e₄(1−κG)", [e1, e4, e1, vm(e4, lin(k4, G))]),
 ("sites 1,2: e₁(1−κ₁Gᵀ)⊗e₄(1−κ₂G)", [vm(e1, lin(k1, GT)), vm(e4, lin(k2, G)), e1, e4]),
 ("pair 1 common g = 1−κGᵀ (A: v g, F: u g^-T)", [vm(e1, lin(k1, GT)), vm(e4, inv(lin(k1, G))), e1, e4]),
 ("pair 2 common g", [e1, e4, vm(e1, lin(k3, GT)), vm(e4, inv(lin(k3, G)))]),
 ("pair 1 common g₁, pair 2 common g₂", [vm(e1, lin(k1, GT)), vm(e4, inv(lin(k1, G))), vm(e1, lin(k3, GT)), vm(e4, inv(lin(k3, G)))]),
 ("global common g (both pairs, same κ)", [vm(e1, lin(k1, GT)), vm(e4, inv(lin(k1, G))), vm(e1, lin(k1, GT)), vm(e4, inv(lin(k1, G)))]),
]
for dp in range(3):
    for dq in range(3):
        if dp + dq <= 2:
            fam.append((f"flag (deg p, deg q) = ({dp},{dq}), block 1 only", flag(dp, dq, 0) + [e1, e4]))
            fam.append((f"flag ({dp},{dq}), both blocks (independent coefficients)", flag(dp, dq, 0) + flag(dp, dq, Fr(1, 3))))
print('| bra (companion frame) | v·u per pair | # c_x ≠ 0 | rank c_(12),(34) | site ranks | rank (13)\|(24) |')
print('|---|---|---|---|---|---|')
for name, vs in fam:
    vs = [X.vec(v) for v in vs]
    vu = [int((vs[0] * vs[1] % p).sum() % p), int((vs[2] * vs[3] % p).sum() % p)]
    vu = ['0' if x == 0 else '≠0' for x in vu]
    r = X.ranks(X.coeffs(X.bra(vs)))
    print(f"| {name} | {vu[0]}, {vu[1]} | {r['nnz']} | {r['pair']} | {r['sites']} | {r['r13_24']} |")
