"""Q-dictionary of the pair-point eigenvalues (pair 1) at the ABJM point, polynomial gauge q_i (qf.py conventions)."""
import numpy as np, itertools
import qf

def make_qI(ch, qs, z, I):
    """reduced Casoratian q_I = P_I / F3,  P_I(u) = det[z_i^{2-s} q_i(u - s h)]_{i in I, s=0,1,2}, F3 = prod_F (u-thF)(u-thF-h)."""
    h = ch.h
    thFs = [t for ty, t in zip(ch.types, ch.th) if ty == 'F']
    def P(u):
        return np.linalg.det(np.array([[z[i] ** (2 - s) * qf.polyval(qs[i], u - s * h) for s in range(3)] for i in I]))
    def F3(u):
        return np.prod([(u - t) * (u - t - h) for t in thFs])
    deg = sum(len(qs[i]) - 1 for i in I) + 0
    pts = np.linspace(-2.7, 2.9, 40) + 0.031
    vals = np.array([P(u) / F3(u) for u in pts])
    # fit polynomial of degree <= deg (P has degree <= sum deg q_i); F3 has degree 2*Lp
    dfit = max(deg - 2 * ch.Lp, 0)
    A = np.vander(pts, dfit + 1, increasing=True)
    c, *_ = np.linalg.lstsq(A, vals, rcond=None)
    resid = np.linalg.norm(A @ c - vals) / max(np.linalg.norm(vals), 1e-300)
    return c, resid

def formulas(ch, qs, z, p=1, I=(0, 1, 2)):
    """Q-formulas for (tau_F, x = tau1(u0), tau2(u0), tau3^(u0)) of pair p (here p = 1 only)."""
    h = ch.h; thF = ch.th[2 * p - 1]; thA = ch.th[2 * p - 2]; u0 = thF + h
    V4 = lambda u: qf.V4(ch, u)
    rest_V4 = lambda u: np.prod([(u - t) if ty == 'A' else (u - t - h) for n, (ty, t) in enumerate(zip(ch.types, ch.th)) if n // 2 != p - 1])
    rest_F3 = lambda u: np.prod([(u - t) * (u - t - h) for n, (ty, t) in enumerate(zip(ch.types, ch.th)) if ty == 'F' and n // 2 != p - 1])
    q = lambda i, u: qf.polyval(qs[i], u)
    tauF = [z[j] * V4(thF) * q(j, thF) / q(j, thF - h) for j in range(4)]
    qI, res = make_qI(ch, qs, z, I)
    ZI = np.prod([z[i] for i in I])
    qIv = lambda u: qf.polyval(qI, u)
    def d(cols):
        return np.linalg.det(np.array([[z[i] ** (-j) * q(i, thA - j * h) for j in cols] for i in I]))
    pref = rest_V4(u0) / rest_F3(thF) * ZI ** 4 / qIv(thF)
    x = pref * d([1, 3, 4])
    t2 = pref * V4(thF) * d([1, 2, 4])
    t3h = h * rest_V4(u0) * V4(thF) * V4(thF - h) * ZI * rest_F3(u0) / rest_F3(thF) * qIv(u0) / qIv(thF)
    return tauF, x, t2, t3h, res
