"""Common python port of the notebook definitions (same H, same vectors) for pasting results."""
import numpy as np
from core import *

def H_diag(chd):
    return chd.t(1, Fr(37, 100)) - 0.7 * chd.t(2, Fr(-113, 100)) + 0.3 * chd.t(3, Fr(3, 2)) + 1.1 * chd.t(1, Fr(-113, 100))

def weights():
    out = []
    for idx in wordsIdx:
        w = [0] * 4
        for s, i in enumerate(idx):
            w[i] += 1 if types[s] == 'A' else -1
        out.append(tuple(w))
    return out
WT = weights()

def states_diag(chd):
    H = H_diag(chd)
    secs = {}
    for k, w in enumerate(WT): secs.setdefault(w, []).append(k)
    V = np.zeros((256, 256), dtype=complex); gaps = []; col = 0
    for w, idx in secs.items():
        ev, U = np.linalg.eig(H[np.ix_(idx, idx)])
        if len(idx) > 1:
            d = np.abs(ev[:, None] - ev[None, :]) + np.eye(len(idx)) * 1e300
            gaps.append(d.min())
        for k in range(len(idx)):
            V[idx, col] = U[:, k]; col += 1
    return V, min(gaps)

def Rmat(z):
    Vf = np.array([[float(x) for x in row] for row in Vand(z)]); VmT = np.linalg.inv(Vf).T
    return np.kron(np.kron(np.kron(VmT, Vf), VmT), Vf)

def tau(op, V):
    TV = op @ V
    return np.einsum('ij,ij->j', V.conj(), TV) / np.einsum('ij,ij->j', V.conj(), V)

def wvals(ch, V):
    th, h = ch.th, ch.h
    n = V.shape[1]
    A1 = np.array([np.ones(n)] + [tau(ch.t(a, th[0] - h), V) for a in (1, 2, 3)])
    A3 = np.array([np.ones(n)] + [tau(ch.t(a, th[2] - h), V) for a in (1, 2, 3)])
    r2 = tau(ch.t(1, th[1]), V); r4 = tau(ch.t(1, th[3]), V)
    return A1, r2, A3, r4

def wprod(Wv):
    A1, r2, A3, r4 = Wv
    return np.array([A1[x[0]] * r2 ** x[1] * A3[x[2]] * r4 ** x[3] for x in wordsIdx])
