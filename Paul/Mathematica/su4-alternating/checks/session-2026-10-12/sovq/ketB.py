import itertools, pickle
import numpy as np
exec(open('bop.py').read().split("W = pickle.load")[0])
ge = pickle.load(open('gram_exact.pkl', 'rb')); Y = np.array([[float(x) for x in r] for r in ge['Y']])
B1 = Bop(0.37); B2 = Bop(-1.13)
def perm_index(rev_idx, rev_sites):
    P = np.zeros((256, 256))
    for n in range(256):
        d = [(n >> (2 * (3 - s))) & 3 for s in range(4)]
        if rev_idx: d = [3 - x for x in d]
        if rev_sites: d = d[::-1]
        m = ((d[0] * 4 + d[1]) * 4 + d[2]) * 4 + d[3]
        P[m, n] = 1
    return P
for ri in (False, True):
    for rs in (False, True):
        P = perm_index(ri, rs)
        cnt_bra = cnt_ket = 0
        for j in range(256):
            y = P @ Y[:, j]
            for B in (None,):
                pass
            # as bra: y^T B ∝ y^T ; as ket: B y ∝ y
            def eig(v, op):
                w = op(v); p = np.argmax(np.abs(v)); lam = w[p] / v[p]
                return np.abs(w - lam * v).max() / max(np.abs(w).max(), 1e-300) < 1e-8
            if eig(y, lambda v: v @ B1) and eig(y, lambda v: v @ B2): cnt_bra += 1
            if eig(y, lambda v: B1 @ v) and eig(y, lambda v: B2 @ v): cnt_ket += 1
        print('index reversal', ri, 'site reversal', rs, ': # ket words that are B-eigen as bras', cnt_bra, ', as kets (B|y> ~ |y>)', cnt_ket)
