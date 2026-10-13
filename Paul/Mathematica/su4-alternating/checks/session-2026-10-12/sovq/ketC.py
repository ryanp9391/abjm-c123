import itertools, pickle
import numpy as np
exec(open('bop.py').read().split("W = pickle.load")[0])
ge = pickle.load(open('gram_exact.pkl', 'rb')); Y = np.array([[float(x) for x in r] for r in ge['Y']])
def qmin_shift(rows, cols, v, sgn):
    a = len(cols); tot = np.zeros((256, 256))
    for p in itertools.permutations(range(a)):
        M = np.identity(256)
        for k in range(a):
            M = M @ T(v + sgn * k * hf)[rows[p[k]]][cols[k]]
        tot += sgn_(p) * M
    return tot
def sgn_(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s
def Cop(u, colsgn, argsgn, rowset):
    tot = np.zeros((256, 256))
    for j1 in rowset:
        for J2 in itertools.combinations(rowset, 2):
            A = T(u)[j1][3]
            Bm = qmin_shift(list(J2), [j1, 3], u + argsgn * hf, colsgn)
            C = qmin_shift(list(rowset), [J2[0], J2[1], 3], u + 2 * argsgn * hf, colsgn)
            for order in ('CBA',):
                tot += C @ Bm @ A
    return tot
def count(op1, op2, kets=True):
    c = 0
    for j in range(256):
        y = Y[:, j]
        ok = True
        for op in (op1, op2):
            w = op @ y if kets else y @ op
            p = np.argmax(np.abs(y)); lam = w[p] / y[p]
            ok &= np.abs(w - lam * y).max() / max(np.abs(w).max(), 1e-300) < 1e-8
        c += ok
    return c
for rowset in ([0, 1, 2], [1, 2, 3]):
    for colsgn in (-1, 1):
        for argsgn in (-1, 1):
            C1 = Cop(0.37, colsgn, argsgn, rowset); C2 = Cop(-1.13, colsgn, argsgn, rowset)
            comm = np.abs(C1 @ C2 - C2 @ C1).max() / np.abs(C1).max() / np.abs(C2).max()
            print('rows', rowset, 'col shifts', colsgn, 'arg shifts', argsgn, '[C,C]=%.1e' % comm, '# ket words C-eigen:', count(C1, C2), ' as bras:', count(C1, C2, False))
