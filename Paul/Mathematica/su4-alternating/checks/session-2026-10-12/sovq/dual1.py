import itertools, pickle
import numpy as np
from sov_l4 import *
hf = 0.5; thf = tuple(float(x) for x in TH0); zf = tuple(float(x) for x in Z0)
ch = Chain(h=hf, th=thf, z=zf, exact=False)
Wd = pickle.load(open('words_exact.pkl', 'rb'))['words']
W = np.array([[float(x) for x in r] for r in Wd])
print('cond(W)', np.linalg.cond(W))
def ketwords(vac, ptsA, ptsF, modeA, modeF):
    """|y> = OA1 OF2 OA3 OF4 |vac>, labels (m1,n2,m3,n4) in 0..3, order irrelevant (commuting)"""
    def op(site, lab):
        th = thf[site]
        if site in (0, 2):
            if modeA == 'pow1': return np.linalg.matrix_power(ch.t(1, th + ptsA), lab)
            if modeA == 'pow3': return np.linalg.matrix_power(ch.t(3, th + ptsA), lab)
            if modeA == 'ta': return ch.t(lab, th + ptsA) if lab else np.identity(256)
            if modeA == 'tastar': return ch.t(lab, th + ptsA + lab * hf) if lab else np.identity(256)
        else:
            if modeF == 'pow1': return np.linalg.matrix_power(ch.t(1, th + ptsF), lab)
            if modeF == 'pow3': return np.linalg.matrix_power(ch.t(3, th + ptsF), lab)
            if modeF == 'ta': return ch.t(lab, th + ptsF) if lab else np.identity(256)
            if modeF == 'tastar': return ch.t(lab, th + ptsF + lab * hf) if lab else np.identity(256)
    cols = []
    for lab in WORDS:
        v = vac.copy()
        for s in (3, 2, 1, 0):
            v = op(s, lab[s]) @ v
        cols.append(v)
    return np.array(cols).T   # columns = kets
def prodket(sites):
    v = np.array([1.0])
    for i in sites:
        e = np.zeros(4); e[i] = 1; v = np.kron(v, e)
    return v
results = []
for vs in [(0, 3, 0, 3), (3, 0, 3, 0)]:
    vac = prodket(vs)
    for modeA, pA in [('pow1', 0.0), ('pow3', -hf), ('ta', -hf), ('pow1', -hf), ('tastar', 0.0), ('tastar', hf)]:
        for modeF, pF in [('tastar', hf), ('tastar', 0.0), ('ta', 0.0), ('pow1', 0.0), ('ta', hf), ('pow3', 0.0), ('tastar', -hf)]:
            Y = ketwords(vac, pA, pF, modeA, modeF)
            rk = np.linalg.matrix_rank(Y, tol=1e-9 * np.abs(Y).max())
            if rk < 256:
                results.append((vs, modeA, pA, modeF, pF, rk, None, None)); continue
            G = W @ Y
            G = G / np.abs(G).max()
            nnz = int((np.abs(G) > 1e-10).sum())
            M = np.linalg.inv(G)
            nnzM = int((np.abs(M) > 1e-10 * np.abs(M).max()).sum())
            results.append((vs, modeA, pA, modeF, pF, rk, nnz, nnzM))
for r in sorted(results, key=lambda r: (r[6] is None, r[6] or 0)):
    print(r)
