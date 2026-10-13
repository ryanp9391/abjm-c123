import sys, pickle
from fractions import Fraction as Fr
import numpy as np
from cyc import modmat_mul, ModBasis
from modrun import Wex_pair, vec_mod
from lib import Z4
Lb_, th_ = int(sys.argv[1]), int(sys.argv[2]); maxl = int(sys.argv[3])
d = np.load(f'cyc_L{Lb_}_th{th_}_p0.npz'); names, Tc, p = pickle.load(open(f'cyc_L{Lb_}_th{th_}_p0.pkl', 'rb'))
ops, e0, rows, piv = d['ops'].astype(np.float64), d['e0'], d['rows'], d['piv']
r = ops.shape[1]; nF = len(names); L = Lb_
a0 = Fr(1, 3); z = Z4
v, u = Wex_pair(a0)
rest = []
for _ in range(L // 2 - 1):
    vv, uu = Wex_pair(a0); rest += [vv, uu]
def pair_target(x, y):   # x: list of 4 (A site), y: list of 4 (F site)
    return vec_mod([x, y] + rest, p)
def singlet(wA):
    # sum_i wA_i e_i (x) e_i ; as a vector in pair-1 space 16-dim, times rest
    vecs = []
    out = np.zeros(16, dtype=object)
    for i in range(4): out[i * 4 + i] = wA[i]
    return out
def mod_of(Fx): return np.array([int(Fr(x).numerator % p) * pow(int(Fr(x).denominator % p), -1, p) % p for x in Fx], dtype=np.float64)
def kron_fr(a, b): return [x * y for x in a for y in b]
restv = [Fr(1)]
for vs in rest: restv = kron_fr(restv, vs)
def pair_full(p16): return np.mod(np.array([ (Fr(x)*y) for x in p16 for y in restv ], dtype=object).astype(object), 1) if False else None
def full_from_pair(p16):
    allv = kron_fr(list(p16), restv)
    return mod_of(allv)
mul = lambda vv, w: [vv[i] * w[i] for i in range(4)]
Lam = z; Laminv = [1 / x for x in z]
extras = {
 'F-site Lam   : v x uLam': full_from_pair(kron_fr(v, mul(u, Lam))),
 'F-site Lam^2 : v x uLam^2': full_from_pair(kron_fr(v, mul(mul(u, Lam), Lam))),
 'F-site Lam^3 : v x uLam^3': full_from_pair(kron_fr(v, mul(mul(mul(u, Lam), Lam), Lam))),
 'F-site Lam^-1: v x uLam^-1': full_from_pair(kron_fr(v, mul(u, Laminv))),
 'A-site Lam   : vLam x u': full_from_pair(kron_fr(mul(v, Lam), u)),
 'A-site Lam^-1: vLam^-1 x u': full_from_pair(kron_fr(mul(v, Laminv), u)),
 'vLam^-1 x uLam': full_from_pair(kron_fr(mul(v, Laminv), mul(u, Lam))),
 'vLam x uLam^-1': full_from_pair(kron_fr(mul(v, Lam), mul(u, Laminv))),
 'singlet S = sum e_i x e_i': full_from_pair(singlet([Fr(1)] * 4)),
 'rho1(Lam^-1) S': full_from_pair(singlet(Laminv)),
 'rho1(Lam) S': full_from_pair(singlet(Lam)),
 'e1 x e4 (vac pair)': full_from_pair(kron_fr([Fr(1), 0, 0, 0], [0, 0, 0, Fr(1)])),
}
ok = {}
for k, vec in extras.items():
    c = vec[piv]
    back = modmat_mul(c[None, :], rows.astype(np.float64), p)[0]
    ok[k] = np.array_equal(back, np.mod(vec, p))
    extras[k] = c
print('in cyclic space:', ok)
B = ModBasis(r, p); B.add(e0[None, :].astype(np.float64).copy()); new = e0[None, :].astype(np.float64)
def rep(l):
    print(f'layer {l}: rank {B.rank} ', ' '.join(('Y' if (ok[k] and B.contains(extras[k])) else 'n') for k in extras), flush=True)
print('order:', list(extras.keys()))
rep(0)
for l in range(1, maxl + 1):
    cand = np.concatenate([modmat_mul(new, ops[i], p) for i in range(nF)], axis=0)
    added = []
    for c0 in range(0, cand.shape[0], 400):
        ch = cand[c0:c0 + 400]; fl = B.add(ch); added.append(ch[np.array(fl)])
    new = np.concatenate(added, axis=0); rep(l)
    if new.shape[0] == 0: break
