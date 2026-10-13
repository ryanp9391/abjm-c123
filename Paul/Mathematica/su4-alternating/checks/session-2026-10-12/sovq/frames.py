"""companion-frame eigenkets / left eigenbras (mp) from the diagonal-frame block eigenvectors"""
import pickle
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
from sov_l4 import vandermonde, minv, mT, Z0, WORDS
mp.mp.dps = 50
tomp = lambda x: mp.mpf(x.numerator) / x.denominator
V = vandermonde(Z0); VmT = mT(minv(V)); Vinv = minv(V); VT = mT(V)
def site_mats(): return VmT, V
lr = pickle.load(open('lr_vecs.pkl', 'rb'))
def apply_kron(mats, vec, left=False):
    """(m1 x m2 x m3 x m4) . vec  (or vec . (m1 x ...) if left) for 4x4 site matrices, vec length 256 (mp list)"""
    import itertools
    T = {}
    v = list(vec)
    for s in range(4):
        m = mats[s]
        new = [mp.mpf(0)] * 256
        stride = 4 ** (3 - s)
        for idx in range(256):
            if v[idx] == 0: continue
            d = (idx // stride) % 4
            base = idx - d * stride
            for k in range(4):
                c = m[d][k] if left else m[k][d]
                if c != 0:
                    new[base + k * stride] += tomp(c) * v[idx]
        v = new
    return v
def companion_ket(n):
    """|Psi_c> = R |Psi_d>, R = V^{-T} (x) V (x) V^{-T} (x) V"""
    return apply_kron([VmT, V, VmT, V], lr['right'][n])
def companion_bra(n):
    """<PsiL_c| = <PsiL_d| R^{-1}, R^{-1} = V^T (x) V^{-1} (x) ..."""
    return apply_kron([VT, Vinv, VT, Vinv], lr['left'][n], left=True)
