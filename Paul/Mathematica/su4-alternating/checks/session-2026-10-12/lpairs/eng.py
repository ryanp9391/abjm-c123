"""eng.py -- fast transfer matrices t_k(u), k = 1,2,3 (fused, as in pgbase.py), DIAGONAL twist Lambda = diag(z),
acting on BRAS psi -> psi . t_k(u), for the alternating chain A F A F ... (site 1 leftmost, = most significant index).
Two arithmetic modes:
   'p' : exact modulo the prime P (all inputs Fractions/ints, converted mod P); identities checked mod P are exact
         up to a false-positive probability ~ 1/P per entry.
   'c' : complex128 numerics.
Conventions = pgbase.py: L^(k)_s = (u - th~_s) + h sum_ab D_k(e_ab) (x) rho_s(E_ab), rho_A(E_ab) = e_ab, rho_F(E_ab) = -e_ba,
th~_A = th, th~_F = th + (k-1)h; prefactor g_A = prod_{j=1}^{k-1}(u-th-jh), g_F = prod_{j=0}^{k-2}(u-th-jh);
'hat' (reg=True): exactly vanishing scalar factors dropped."""
import numpy as np, itertools
from fractions import Fraction as Fr
from pgbase import DK, subsets

P = 268435399

def modp(x, p=P):
    x = Fr(x)
    return (x.numerator % p) * pow(x.denominator % p, p - 2, p) % p

_S = {}
def struct(k, ty):
    if (k, ty) in _S: return _S[(k, ty)]
    n = len(subsets(k)); S = np.zeros((n, n, 4, 4), dtype=np.int64)
    for a in range(4):
        for b in range(4):
            D = np.array(DK[k][(a, b)], dtype=np.int64)
            R = np.zeros((4, 4), dtype=np.int64)
            if ty == 'A': R[a, b] = 1
            else: R[b, a] = -1
            S += np.einsum('KJ,pq->KJpq', D, R)
    _S[(k, ty)] = S
    return S

def gfactors(ty, k, u, th, h):
    js = range(1, k) if ty == 'A' else range(0, k - 1)
    return [u - th - j * h for j in js]

def iszero(x, mode):
    if mode == 'p': return x == 0
    return abs(x) < 1e-11

class Chain:
    def __init__(self, Lp, th, h, z, mode='c'):
        self.Lp = Lp; self.N = 2 * Lp; self.types = ['A', 'F'] * Lp
        self.th = list(th); self.h = h; self.z = list(z); self.mode = mode
    def cv(self, x):
        """convert scalar to the working field"""
        if self.mode == 'p': return modp(x)
        return complex(x)
    def wedge_diag(self, k):
        out = []
        for I in subsets(k):
            p = 1
            for i in I: p = p * self.z[i]
            out.append(p)
        return out
    def pref(self, k, u, reg):
        r = 1
        for ty, t in zip(self.types, self.th):
            for f in gfactors(ty, k, u, t, self.h):
                if reg and iszero(f, self.mode): continue
                r = r * f
        return r
    def apply_t(self, k, u, psi, reg=False):
        """psi: array (4,)*N or (4,)*N + (B,).  Returns psi . t_k(u) (same shape)."""
        N = self.N; h = self.h; mode = self.mode
        squeeze = (psi.ndim == N)
        if squeeze: psi = psi[..., None]
        n = len(subsets(k))
        dt = np.int64 if mode == 'p' else complex
        Phi = np.zeros((n, n) + psi.shape, dtype=dt)
        for I in range(n): Phi[I, I] = psi
        hh = self.cv(h)
        for s in range(N):
            ty = self.types[s]; t = self.th[s]
            tt = t if ty == 'A' else t + (k - 1) * h
            c0 = self.cv(u - tt)
            S = struct(k, ty)
            L = np.zeros((n, n, 4, 4), dtype=dt)
            for K in range(n):
                for p in range(4): L[K, K, p, p] = c0
            if mode == 'p':
                L = (L + hh * (S % P)) % P
            else:
                L = L + hh * S
            ax = 2 + s
            Phi = np.tensordot(Phi, L, axes=([1, ax], [0, 2]))
            Phi = np.moveaxis(Phi, -2, 1)
            Phi = np.moveaxis(Phi, -1, ax)
            if mode == 'p': Phi %= P
        wg = [self.cv(x) for x in self.wedge_diag(k)]
        out = np.zeros(psi.shape, dtype=dt)
        for I in range(n):
            if mode == 'p': out = (out + wg[I] * Phi[I, I]) % P
            else: out = out + wg[I] * Phi[I, I]
        pf = self.cv(self.pref(k, u, reg))
        out = (pf * out) % P if mode == 'p' else pf * out
        return out[..., 0] if squeeze else out

def prod_state(vecs, mode='c'):
    if mode == 'p':
        vv = [np.array([modp(x) for x in v], dtype=np.int64) for v in vecs]
        out = vv[0]
        for v in vv[1:]: out = np.multiply.outer(out, v) % P
        return out
    out = np.array(vecs[0], dtype=complex)
    for v in vecs[1:]: out = np.multiply.outer(out, np.array(v, dtype=complex))
    return out

def W(a):
    return [[1, 0, 0, a], [-a, 0, 0, 1]]

VAC = [[1, 0, 0, 0], [0, 0, 0, 1]]

def lam(ch, k, u, p):
    """localisation constant of t_k(u) on pair p (1-based): product over pairs n != p (vanishing factors dropped)"""
    h = ch.h; r = 1
    for n in range(ch.Lp):
        if n == p - 1: continue
        tA = ch.th[2 * n]; tF = ch.th[2 * n + 1]
        fac = [u - tA - j * h for j in range(k)] + [u - tF - k * h] + [u - tF - j * h for j in range(k - 1)]
        for f in fac:
            if iszero(f, ch.mode): continue
            r = r * f
    return r

def comb(ch, coeffs, vecs):
    """sum c_i v_i in the working field (coeffs given as Fractions/complex)"""
    if ch.mode == 'p':
        out = np.zeros_like(vecs[0])
        for c, v in zip(coeffs, vecs): out = (out + modp(c) * v) % P
        return out
    out = np.zeros_like(vecs[0], dtype=complex)
    for c, v in zip(coeffs, vecs): out = out + complex(c) * v
    return out
