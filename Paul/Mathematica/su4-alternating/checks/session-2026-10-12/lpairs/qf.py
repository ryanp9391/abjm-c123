"""qf.py -- eigenvalues and twisted Q-functions of the alternating chain (numeric, complex128), any Lp (small).
Baxter (sov_l4 convention): sum_{a=0}^4 (-1)^a tau_a(u) Q(u - a h) = 0, Q_i(u) = z_i^{u/h} Phi(u) q_i(u),
Phi(u) = V4(u) Phi(u-h), V4(u) = prod_A (u - thA) prod_F (u - thF - h):
   sum_a (-1)^a z_i^{-a} tau_a(u) q_i(u - a h) prod_{b=a}^{3} V4(u - b h) = 0   (polynomial form)."""
import numpy as np, itertools
import eng

def full_matrix(ch, k, u, reg=False):
    D = 4 ** ch.N
    E = np.eye(D, dtype=complex).reshape((4,) * ch.N + (D,))
    return ch.apply_t(k, u, E, reg=reg).reshape(D, D)   # M[p,q] = (e_p . t)_q ; ket action M @ psi

def joint_eigen(ch, gen=((1.0, 1, 0.37, False), (0.61, 3, -0.29, True), (0.23, 2, 0.83, False))):
    M = sum(c * full_matrix(ch, k, u, reg) for c, k, u, reg in gen)
    ev, V = np.linalg.eig(M)
    gaps = np.array([np.min(np.abs(np.delete(ev, i) - ev[i])) for i in range(len(ev))])
    return V, ev, gaps

def eigvals_at(ch, V, k, u, reg=False):
    M = full_matrix(ch, k, u, reg)
    MV = M @ V
    tau = np.einsum('ij,ij->j', V.conj(), MV) / np.einsum('ij,ij->j', V.conj(), V)
    res = np.linalg.norm(MV - V * tau[None, :], axis=0) / np.maximum(np.abs(tau), 1e-300)
    return tau, res

def V4(ch, u):
    r = 1.0
    for ty, t in zip(ch.types, ch.th):
        r *= (u - t) if ty == 'A' else (u - t - ch.h)
    return r

def baxter_rows(ch, taus_at, us, zi, deg):
    """rows of the linear system for q's coefficients (q(u) = sum c_m u^m, m = 0..deg) at sample points us.
    taus_at[a][n] = tau_a(us[n]) for one state (a = 0..4)."""
    h = ch.h; rows = []
    for n, u in enumerate(us):
        row = np.zeros(deg + 1, dtype=complex)
        for a in range(5):
            pref = (-1) ** a * zi ** (-a) * taus_at[a][n]
            for b in range(a, 4): pref = pref * V4(ch, u - b * h)
            row += pref * (u - a * h) ** np.arange(deg + 1)
        rows.append(row)
    A = np.array(rows)
    sc = np.abs(A).max(axis=1, keepdims=True); sc[sc == 0] = 1
    return A / sc

def solve_q(ch, taus_at, us, zi, maxdeg):
    """lowest degree monic solution; returns (coeffs low->high, smallest singular value ratio)"""
    best = None
    for deg in range(0, maxdeg + 1):
        A = baxter_rows(ch, taus_at, us, zi, deg)
        u_, s, vh = np.linalg.svd(A)
        rat = s[-1] / s[0]
        if rat < 1e-9:
            c = vh[-1].conj(); c = c / c[-1]
            return c, rat
        best = (deg, rat)
    return None, best

def polyval(c, u):
    return sum(cm * u ** m for m, cm in enumerate(c))
