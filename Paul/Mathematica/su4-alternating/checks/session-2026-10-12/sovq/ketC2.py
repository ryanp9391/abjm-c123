import itertools, pickle
import numpy as np
exec(open('bop.py').read().split("W = pickle.load")[0])
ge = pickle.load(open('gram_exact.pkl', 'rb')); Y = np.array([[float(x) for x in r] for r in ge['Y']])
def Cstar(u, c2, c3):
    tot = np.zeros((256, 256))
    for k1 in (1, 2, 3):
        for K2 in itertools.combinations((1, 2, 3), 2):
            tot += qmin([0, 1, 2], [0, K2[0], K2[1]], u + c3) @ qmin([k - 1 for k in K2], [0, k1], u + c2) @ T(u)[k1 - 1][0]
    return tot
def count(op1, op2, kets=True):
    c = 0
    for j in range(256):
        y = Y[:, j]; ok = True
        for op in (op1, op2):
            w = op @ y if kets else y @ op
            if np.abs(w).max() < 1e-12 * np.abs(op).max(): ok = False; break
            p = np.argmax(np.abs(y)); lam = w[p] / y[p]
            ok &= np.abs(w - lam * y).max() / np.abs(w).max() < 1e-8
        c += ok
    return c
for c2, c3 in [(-hf, -2 * hf), (0, -hf), (hf, 2 * hf), (0, hf), (-hf, 0), (hf, 0), (-2 * hf, -3 * hf), (hf, hf)]:
    C1 = Cstar(0.37, c2, c3); C2 = Cstar(-1.13, c2, c3)
    nrm = np.abs(C1).max()
    comm = np.abs(C1 @ C2 - C2 @ C1).max() / (nrm * np.abs(C2).max())
    print('shifts', (c2, c3), '|C|=%.2e' % nrm, '[C,C]=%.1e' % comm, '# ket words C-eigen:', count(C1, C2), ' bra words?')
