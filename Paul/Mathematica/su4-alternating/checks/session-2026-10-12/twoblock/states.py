"""Twisted eigenstates (float), computed in the diagonal frame weight sector by weight sector,
mapped to the companion frame by R (Psi_comp = R Psi_diag). tau-values at the SoV points."""
import numpy as np
from core import *


def weights():
    W = []
    for idx in wordsIdx:  # same ordering as the product basis (i1,i2,i3,i4), 0-based
        w = [0] * 4
        for s, i in enumerate(idx):
            w[i] += 1 if types[s] == 'A' else -1
        W.append(tuple(w))
    return W


def Rmat(z, mode='float'):
    """R = (V^{-T} on A) x V (F) x V^{-T} x V ;  <x_diag| = <x_comp| R ; Psi_comp = R Psi_diag"""
    V = np.array([[float(x) for x in row] for row in Vand(z)])
    VmT = np.linalg.inv(V).T
    return np.kron(np.kron(np.kron(VmT, V), VmT), V)


def eigenstates_diag(chd, seed=3):
    """chd: Chain with diagonal twist. Returns kets (columns) and eigen-values of a generic combination."""
    rng = np.random.default_rng(seed)
    ops = [chd.t(a, u) for a in (1, 2, 3) for u in (Fr(37, 100), Fr(-113, 100), Fr(3, 2))]
    c = rng.normal(size=len(ops))
    H = sum(ci * o for ci, o in zip(c, ops))
    W = weights()
    sectors = {}
    for k, w in enumerate(W):
        sectors.setdefault(w, []).append(k)
    Vs = np.zeros((256, 256), dtype=complex)
    evs = np.zeros(256, dtype=complex)
    col = 0
    for w, idx in sectors.items():
        Hs = H[np.ix_(idx, idx)]
        ev, U = np.linalg.eig(Hs)
        for k in range(len(idx)):
            Vs[idx, col] = U[:, k]
            evs[col] = ev[k]
            col += 1
    return evs, Vs


def tau(ch, a, u, V):
    if a == 0:
        return np.ones(V.shape[1], dtype=complex)
    TV = ch.t(a, u) @ V
    return np.einsum('ij,ij->j', V.conj(), TV) / np.einsum('ij,ij->j', V.conj(), V)


def wavevals(ch, V):
    th, h = ch.th, ch.h
    WA1 = np.array([tau(ch, a, th[0] - h, V) for a in range(4)])
    WA3 = np.array([tau(ch, a, th[2] - h, V) for a in range(4)])
    t2 = tau(ch, 1, th[1], V)
    t4 = tau(ch, 1, th[3], V)
    WF2 = np.array([t2 ** k for k in range(4)])
    WF4 = np.array([t4 ** k for k in range(4)])
    return WA1, WF2, WA3, WF4


def Wprod(Ws):
    WA1, WF2, WA3, WF4 = Ws
    return np.einsum('as,bs,cs,ds->abcds', WA1, WF2, WA3, WF4).reshape(256, -1)
