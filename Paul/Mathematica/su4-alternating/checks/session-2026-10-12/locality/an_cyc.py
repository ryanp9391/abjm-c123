import sys, pickle, itertools, time
import numpy as np
from cyc import modmat_mul, mod_inv_matrix, ModBasis
Lb_, th_ = int(sys.argv[1]), int(sys.argv[2])
d = np.load(f'cyc_L{Lb_}_th{th_}_p0.npz'); names, Tc, p = pickle.load(open(f'cyc_L{Lb_}_th{th_}_p0.pkl', 'rb'))
ops, e0, Lb = d['ops'], d['e0'], d['Lb']
r = ops.shape[1]; nF = len(names)
print('L', Lb_, 'theta pt', th_, 'cyclic dim', r, 'loc (basis rows)', Lb.shape[0], flush=True)
T1 = [k for k in Tc if k.startswith('T1')]; T2 = [k for k in Tc if k.startswith('T2')]; T3 = [k for k in Tc if k.startswith('T3')]
Aidx = [i for i, n in enumerate(names) if n[1] == 0]; Fidx = [i for i, n in enumerate(names) if n[1] == 1]
def report(label, rows):
    B = ModBasis(r, p)
    for c0 in range(0, rows.shape[0], 500): B.add(rows[c0:c0 + 500].copy())
    rk = B.rank
    # loc dim = rank(span)+rank(Lb)-rank(span+Lb)
    B2 = ModBasis(r, p); B2.add(Lb.copy()); rl = B2.rank
    B3 = ModBasis(r, p)
    allr = np.concatenate([rows, Lb], axis=0)
    for c0 in range(0, allr.shape[0], 500): B3.add(allr[c0:c0 + 500].copy())
    ld = rk + rl - B3.rank
    t = lambda ks: ''.join('Y' if B.contains(Tc[k]) else 'n' for k in ks)
    print(f'{label:46s} nwords {rows.shape[0]:6d} rank {rk:4d} loc {ld} T1:{t(T1)} T2:{t(T2)} T3:{t(T3)}', flush=True)
    return B
# word generation
def mul(v, i): return modmat_mul(v, ops[i].astype(np.float64), p)
e0r = e0[None, :].astype(np.float64)
W1 = [e0r] + [mul(e0r, i) for i in range(nF)]
single = {i: mul(e0r, i) for i in range(nF)}
pair = {}
for i in range(nF):
    row_i = single[i]
    block = mul(row_i, 0) if False else None
for i in range(nF):
    for j in range(i, nF):
        pair[(i, j)] = mul(single[i], j)
rows1 = np.concatenate([e0r] + [single[i] for i in range(nF)], axis=0)
report('<=1 factor (all types)', rows1)
rows2 = np.concatenate([rows1] + [pair[k] for k in pair], axis=0)
report('<=2 factors (all types)', rows2)
# triples AFF, AAF, + FFF, AAA
def triples(types_pattern):
    out = []
    for (i, j), v in pair.items():
        for k in range(j, nF):
            pass
    return out
trip = []
labs3 = []
for i in range(nF):
    for j in range(i, nF):
        for k in range(j, nF):
            labs3.append((i, j, k))
# vectorised: group by (i,j) -> apply all k>=j at once
t0 = time.time()
rows3 = []
cls = []
for (i, j), v in pair.items():
    ks = list(range(j, nF))
    V = np.concatenate([mul(v, k) for k in ks], axis=0)
    rows3.append(V)
    for k in ks:
        ty = sorted(['A' if names[x][1] == 0 else 'F' for x in (i, j, k)])
        cls.append(''.join(ty))
rows3 = np.concatenate(rows3, axis=0); cls = np.array(cls)
print('triples built', rows3.shape, time.time() - t0, flush=True)
selAFF = np.isin(cls, ['AFF', 'AAF'])
report('<=2 + triples AFF/AAF (task classes)', np.concatenate([rows2, rows3[selAFF]], axis=0))
report('<=2 + triples FFF', np.concatenate([rows2, rows3[cls == 'FFF']], axis=0))
report('<=2 + triples AAA', np.concatenate([rows2, rows3[cls == 'AAA']], axis=0))
report('<=3 (all types)', np.concatenate([rows2, rows3], axis=0))
# W4 class: (1, t1(th1), t1(th1)^-1) x t1(th2+m1 h) t2(th2+m2 h) t3(th2+m3 h)
idx = {n: i for i, n in enumerate(names)}
t1th1 = idx[(1, 0, 0)]
inv_t1th1 = mod_inv_matrix(ops[t1th1].astype(np.float64), p)
W4rows = {'1': [], 't1(th1)': [], 't1(th1)^-1': []}
for m1 in range(-4, 5):
    for m2 in range(-4, 5):
        for m3 in range(-4, 5):
            if (1, 1, m1) in idx and (2, 1, m2) in idx and (3, 1, m3) in idx:
                v = mul(mul(mul(e0r, idx[(1, 1, m1)]), idx[(2, 1, m2)]), idx[(3, 1, m3)])
                W4rows['1'].append(v)
                W4rows['t1(th1)'].append(mul(v, t1th1))
                W4rows['t1(th1)^-1'].append(modmat_mul(v, inv_t1th1, p))
W4 = {k: np.concatenate(v, axis=0) for k, v in W4rows.items()}
base3 = np.concatenate([rows2, rows3[selAFF]], axis=0)
report('task classes + FFF t1t2t3 (x1)', np.concatenate([base3, W4['1']], axis=0))
report('task classes + (x t1(th1))', np.concatenate([base3, W4['1'], W4['t1(th1)']], axis=0))
report('task classes + (x t1(th1)^-1)', np.concatenate([base3, W4['1'], W4['t1(th1)^-1']], axis=0))
report('task classes + all W4', np.concatenate([base3, W4['1'], W4['t1(th1)'], W4['t1(th1)^-1']], axis=0))
# inverse words: named inverses and all invertible inverses
def inv_rows(idxs, base_rows):
    out = []
    for i in idxs:
        try: M = mod_inv_matrix(ops[i].astype(np.float64), p)
        except ZeroDivisionError: continue
        out.append(modmat_mul(base_rows, M, p))
    return np.concatenate(out, axis=0)
t1th2 = idx[(1, 1, 0)]
base1 = rows1; base2 = rows2
report('<=2 positive + {t1(th1)^-1, t1(th2)^-1} x (1,singles)', np.concatenate([rows2, inv_rows([t1th1, t1th2], rows1)], axis=0))
report('<=2 positive + {t1(th1)^-1, t1(th2)^-1} x (<=2)', np.concatenate([rows2, inv_rows([t1th1, t1th2], rows2)], axis=0))
rows_inv_all = inv_rows(range(nF), rows1)
report('<=2 positive + all inverses x (1,singles)', np.concatenate([rows2, rows_inv_all], axis=0))
report('<=3 positive + {t1(th1)^-1,t1(th2)^-1} x (<=2)', np.concatenate([rows2, rows3, inv_rows([t1th1, t1th2], rows2)], axis=0))
report('<=3 all + all inverses x (1,singles)', np.concatenate([rows2, rows3, rows_inv_all], axis=0))
