import sys, pickle, itertools, time
import numpy as np
from cyc import modmat_mul, mod_inv_matrix, ModBasis
Lb_, th_ = int(sys.argv[1]), int(sys.argv[2])
d = np.load(f'cyc_L{Lb_}_th{th_}_p0.npz'); names, Tc, p = pickle.load(open(f'cyc_L{Lb_}_th{th_}_p0.pkl', 'rb'))
ops, e0, Lb = d['ops'].astype(np.float64), d['e0'].astype(np.float64), d['Lb']
r = ops.shape[1]; nF = len(names)
t0 = time.time()
inv = []
for i in range(nF):
    try: inv.append(mod_inv_matrix(ops[i], p))
    except ZeroDivisionError: inv.append(None)
print('inverses built', sum(m is not None for m in inv), 'of', nF, time.time() - t0, flush=True)
np.savez_compressed(f'inv_L{Lb_}_th{th_}.npz', inv=np.stack([m if m is not None else np.zeros((r, r)) for m in inv]))
T1 = [k for k in Tc if k.startswith('T1')]; T2 = [k for k in Tc if k.startswith('T2')]; T3 = [k for k in Tc if k.startswith('T3')]
def rep(label, rows):
    B = ModBasis(r, p)
    for c0 in range(0, rows.shape[0], 500): B.add(rows[c0:c0 + 500].copy())
    rk = B.rank
    B2 = ModBasis(r, p); B2.add(Lb.copy()); rl = B2.rank
    B3 = ModBasis(r, p); allr = np.concatenate([rows, Lb], axis=0)
    for c0 in range(0, allr.shape[0], 500): B3.add(allr[c0:c0 + 500].copy())
    ld = rk + rl - B3.rank
    t = lambda ks: ''.join('Y' if B.contains(Tc[k]) else 'n' for k in ks)
    print(f'{label:62s} nwords {rows.shape[0]:6d} rank {rk:4d} loc {ld} T1:{t(T1)} T2:{t(T2)} T3:{t(T3)}', flush=True)
e0r = e0[None, :]
mul = lambda v, M: modmat_mul(v, M, p)
single = np.concatenate([mul(e0r, ops[i]) for i in range(nF)], axis=0)
invs = np.concatenate([mul(e0r, inv[i]) for i in range(nF) if inv[i] is not None], axis=0)
okidx = [i for i in range(nF) if inv[i] is not None]
rep('1 + positive singles', np.concatenate([e0r, single]))
rep('1 + inverse singles', np.concatenate([e0r, invs]))
rep('1 + positive singles + inverse singles', np.concatenate([e0r, single, invs]))
# inverse of single f times single g
ig = np.concatenate([np.concatenate([mul(single[g:g + 1], inv[i]) for g in range(nF)], axis=0) for i in okidx], axis=0)
rep('1 + singles + inv singles + (inv f)(single g)', np.concatenate([e0r, single, invs, ig]))
# by type of the inverse factor
def fam(pred): return [i for i in okidx if pred(names[i])]
for label, idxs in [('inverse t1 only', fam(lambda n: n[0] == 1)), ('inverse t2 only', fam(lambda n: n[0] == 2)), ('inverse t3 only', fam(lambda n: n[0] == 3)),
                    ('inverse A-site (s=1) only', fam(lambda n: n[1] == 0)), ('inverse F-site (s=2) only', fam(lambda n: n[1] == 1)),
                    ('inverse m=0 only', fam(lambda n: n[2] == 0)), ('inverse t1(th1), t1(th2) (task)', [i for i in okidx if names[i] in [(1, 0, 0), (1, 1, 0)]])]:
    rows = np.concatenate([e0r, single] + [mul(e0r, inv[i]) for i in idxs] + [np.concatenate([mul(single[g:g + 1], inv[i]) for g in range(nF)], axis=0) for i in idxs], axis=0)
    rep(f'1 + singles + {label} x (1,single)  [{len(idxs)} inverses]', rows)
