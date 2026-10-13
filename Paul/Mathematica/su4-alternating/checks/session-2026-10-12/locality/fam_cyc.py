import sys, pickle
import numpy as np
from cyc import modmat_mul, ModBasis
Lb_, th_, maxl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
d = np.load(f'cyc_L{Lb_}_th{th_}_p0.npz'); names, Tc, p = pickle.load(open(f'cyc_L{Lb_}_th{th_}_p0.pkl', 'rb'))
ops, e0 = d['ops'].astype(np.float64), d['e0'].astype(np.float64)
r = ops.shape[1]
T12 = [k for k in Tc if k.startswith('T1') or k.startswith('T2')]
def layered(fidx):
    B = ModBasis(r, p); B.add(e0[None, :].copy()); new = e0[None, :]
    out = []
    for l in range(1, maxl + 1):
        cand = np.concatenate([modmat_mul(new, ops[i], p) for i in fidx], axis=0)
        added = []
        for c0 in range(0, cand.shape[0], 400):
            ch = cand[c0:c0 + 400]; fl = B.add(ch); added.append(ch[np.array(fl)])
        new = np.concatenate(added, axis=0)
        out.append((l, B.rank, all(B.contains(Tc[k]) for k in T12)))
        if new.shape[0] == 0: break
    return out
idx = lambda pred: [i for i, n in enumerate(names) if pred(n)]
fam = {
 'all 48': idx(lambda n: True),
 'A-site only': idx(lambda n: n[1] == 0),
 'F-site only': idx(lambda n: n[1] == 1),
 't1 only': idx(lambda n: n[0] == 1), 't2 only': idx(lambda n: n[0] == 2), 't3 only': idx(lambda n: n[0] == 3),
 't1,t3': idx(lambda n: n[0] in (1, 3)), 't1,t2': idx(lambda n: n[0] in (1, 2)), 't2,t3': idx(lambda n: n[0] in (2, 3)),
 '|m|<=2': idx(lambda n: abs(n[2]) <= 2), '|m|<=1': idx(lambda n: abs(n[2]) <= 1), 'm in -2..0': idx(lambda n: -2 <= n[2] <= 0),
 'm in 0..2': idx(lambda n: 0 <= n[2] <= 2), 'm in -1..0': idx(lambda n: -1 <= n[2] <= 0),
 'm=0 only': idx(lambda n: n[2] == 0),
}
for k, fi in fam.items():
    res = layered(fi)
    first = next((l for l, rk, ok in res if ok), None)
    print(f'{k:14s} nfactors {len(fi):2d}  ranks by layer {[rk for l, rk, ok in res]}  first layer with T1,T2: {first}', flush=True)
