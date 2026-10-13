"""lp.py -- pair words K^(p) on the alternating chain at the ABJM point (helpers shared by the task scripts)."""
import numpy as np
from fractions import Fraction as Fr
import eng
from eng import P, modp

def chain_abjm(Lp, thF, h, z, mode):
    th = []
    for f in thF: th += [f + 2 * h, f]
    return eng.Chain(Lp, th, h, z, mode)

def Xs(ch, psi, p):
    """the four non-trivial pair-p words applied to psi: X1, X2, X3, X3X1  (pair p 1-based)"""
    h = ch.h; uF = ch.th[2 * p - 1]; u2 = uF + h
    x1 = ch.apply_t(1, uF, psi)
    x2 = ch.apply_t(2, u2, psi)
    x3 = ch.apply_t(3, u2, psi, reg=True)
    x31 = ch.apply_t(1, uF, x3)
    return x1, x2, x3, x31

def Kcoef(ch, p, kap, y):
    """coefficients of (1, X1, X2, X3, X3X1) in K^(p)_kappa with lambda's of pair p (product over the other pairs)"""
    h = ch.h; uF = ch.th[2 * p - 1]; u2 = uF + h; z = ch.z
    l1, l2, l3 = eng.lam(ch, 1, uF, p), eng.lam(ch, 2, u2, p), eng.lam(ch, 3, u2, p)
    s = 1 - kap * (z[0] + z[3])
    return [s * s, kap * s / (h**2 * l1), kap * s / (2 * h**4 * y * l2), kap * s / (6 * h**5 * y**2 * l3),
            kap**2 / (24 * h**7 * y**2 * l1 * l3)]

def applyK(ch, psi, p, kap, y):
    x = Xs(ch, psi, p)
    return eng.comb(ch, Kcoef(ch, p, kap, y), [psi] + list(x))

def Wstate(ch, alist):
    vecs = []
    for a in alist: vecs += eng.W(a)
    return eng.prod_state(vecs, ch.mode)

def eq(ch, A, B):
    if ch.mode == 'p': return bool(np.array_equal(A % P, B % P))
    return float(np.abs(A - B).max() / max(np.abs(B).max(), 1e-300))
