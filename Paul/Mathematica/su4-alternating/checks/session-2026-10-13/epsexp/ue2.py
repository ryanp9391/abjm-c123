"""ue2.py -- extensions of ue_base.py for the epsilon-expansion study (session 2026-10-13).
 * apply_t1G: t_1^G(u) with an ARBITRARY twist matrix G (bras psi -> psi.t1^G), used for G_a = Lambda +- a (z1-z4) E41 and T_14(u) (= t1 with G = E41)
 * odd-bridge chain (NEXT.md 2j): closed constants d_k, bras P_k, word eigenvalues
 * helpers: <vac| J^k, J = global rho(E14), weights."""
import numpy as np, itertools, sys, os
from ue_base import *

def apply_t1G(ch, u, psi, G):
    """psi . t_1^G(u),  t_1^G = tr_aux[L_1 ... L_N G]   (aux = C^4; G arbitrary 4x4)"""
    N = ch.N; hh = ch.cv(ch.h)
    squeeze = (psi.ndim == N)
    if squeeze: psi = psi[..., None]
    dt = complex
    Phi = np.zeros((4, 4) + psi.shape, dtype=dt)
    for I in range(4): Phi[I, I] = psi
    for s in range(N):
        ty = ch.types[s]; t = ch.th[s]
        c0 = u - t
        S = eng.struct(1, ty)
        L = np.zeros((4, 4, 4, 4), dtype=dt)
        for K in range(4):
            for p in range(4): L[K, K, p, p] = c0
        L = L + hh * S
        ax = 2 + s
        Phi = np.tensordot(Phi, L, axes=([1, ax], [0, 2]))
        Phi = np.moveaxis(Phi, -2, 1)
        Phi = np.moveaxis(Phi, -1, ax)
    out = np.zeros(psi.shape, dtype=dt)
    for I in range(4):
        for J in range(4):
            if G[J, I] != 0: out = out + G[J, I] * Phi[I, J]
    return out[..., 0] if squeeze else out

def vac_bra(S):
    return S.vac()

def e41():
    G = np.zeros((4, 4), dtype=complex); G[3, 0] = 1; return G
def e14():
    G = np.zeros((4, 4), dtype=complex); G[0, 3] = 1; return G

def tau_ket(ch, psi, k, u, reg=False):
    """eigenvalue of t_k(u) on the ket psi, from one row (basis bra) -- psi must be an eigenvector"""
    p = int(np.argmax(np.abs(psi)))
    N = ch.N
    e = np.zeros(4 ** N, dtype=complex); e[p] = 1
    row = ch.apply_t(k, u, e.reshape((4,) * N), reg=reg).reshape(-1)
    return (row @ psi) / psi[p]

# ---- odd-bridge chain constants (ob.py conventions, theta_k k = 1..N, types A,F,A,F,...)
def d_closed(k, th, Lp, h=0.5):
    N = 2 * Lp; ty = ['A', 'F'] * Lp
    if k % 2 == 1:
        p = h
        for s in range(1, N + 1):
            if s == k: continue
            p *= (th[k - 1] - th[s - 1]) if ty[s - 1] == 'F' else (th[k - 1] - th[s - 1] + h)
        return 1 / p
    p = -h
    for s in range(1, N + 1):
        if s == k: continue
        p *= (th[k - 1] - th[s - 1]) if ty[s - 1] == 'A' else (th[k - 1] - th[s - 1] - h)
    return p
def c_closed(k, th, Lp):
    r = 1.0
    for j in range(1, k + 1): r *= d_closed(j, th, Lp)
    return r

def Wv(a): return [np.array([1, 0, 0, a], dtype=complex), np.array([-a, 0, 0, 1], dtype=complex)]
def kron_state(sv):
    out = sv[0]
    for v in sv[1:]: out = np.multiply.outer(out, v)
    return out.reshape(-1)

def chain_bra(Lp, k, a, z, eps, cc):
    """bra of the chain word with last free letter at site k (1-based), sites j<k torus rotated (A: Lambda^-1, F: Lambda), site k free,
    sites > k = W_a (ob.py site_vectors)."""
    Lam = np.diag(z); LamI = np.diag(1 / np.array(z)); I4 = np.eye(4)
    sv = Wv(a) * Lp
    for j in range(k - 1): sv[j] = sv[j] @ (LamI if j % 2 == 0 else Lam)
    if k % 2 == 0: sv[k - 1] = sv[k - 1] @ ((cc * eps - 1) * I4 + Lam)
    else:
        kap = 1 - cc * eps; sv[k - 1] = sv[k - 1] @ (I4 - kap * LamI)
    return kron_state(sv)

def chain_word_eig(Lp, k, th, taus1, cc, eps):
    """eigenvalue of the word prod_{j<k} t1(th_j)^{(-1)^j-ish}(al + be t1(th_k)^{+-1}); taus1[j-1] = tau1(th_j)"""
    if k % 2 == 0:
        al, be = (cc * eps - 1) / c_closed(k - 1, th, Lp), 1 / c_closed(k, th, Lp)
    else:
        kap = 1 - cc * eps; al, be = 1 / c_closed(k - 1, th, Lp), -kap / c_closed(k, th, Lp)
    pre = 1.0
    for j in range(1, k):
        pre *= (1 / taus1[j - 1]) if j % 2 == 1 else taus1[j - 1]
    last = taus1[k - 1]
    return pre * (al + be * (last if k % 2 == 0 else 1 / last))
