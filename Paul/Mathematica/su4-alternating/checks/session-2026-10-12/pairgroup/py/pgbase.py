"""Base module: su(4) alternating chain A F A F ..., fused transfer matrices t_k (k=0..4),
bra (row-vector) conventions psi -> psi.M.  Works with float/complex numpy arrays or object arrays of Fractions."""
import numpy as np, itertools
from fractions import Fraction as Fr

def subsets(k): return list(itertools.combinations(range(4), k))

def minor(M, I, J):
    # determinant of submatrix rows I cols J (generic dtype, small)
    k = len(I)
    if k == 0: return 1
    tot = 0
    for p in itertools.permutations(range(k)):
        s = perm_sign(p); term = s
        for a in range(k): term = term * M[I[a], J[p[a]]]
        tot = tot + term
    return tot

def perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

def wedge(M, k, dtype=object):
    S = subsets(k)
    W = np.empty((len(S), len(S)), dtype=dtype)
    for a, I in enumerate(S):
        for b, J in enumerate(S):
            W[a, b] = minor(M, I, J)
    return W

def eab(a, b, dtype=object, one=1):
    m = np.zeros((4, 4), dtype=dtype); m[a, b] = one; return m

def Dk(k, a, b, dtype=object):
    """derivation action of e_ab on wedge^k (matrix (I,J) = d/de minor_{IJ}(1+e e_ab))"""
    S = subsets(k); n = len(S)
    D = np.zeros((n, n), dtype=dtype)
    if dtype == object: D[:] = 0
    # (1+e e_ab) minors to first order: entry (I,J) of wedge^k(1+eX) = det((1+eX)[I,J])
    # derivative = sum over positions: replace one row... compute directly
    for x, I in enumerate(S):
        for y, J in enumerate(S):
            # d/de det(A + e X)[I,J] at e=0 with A = identity: sum_{r,c} cofactor(Id[I,J])_{rc} X[I_r,J_c]
            # do numerically-exact: det of (Id + e X) restricted, linear coefficient via 2-point formula with integers
            val = 0
            k_ = len(I)
            for r in range(k_):
                for c in range(k_):
                    if I[r] == a and J[c] == b:
                        # cofactor of Id[I,J] at (r,c)
                        rows = [I[t] for t in range(k_) if t != r]; cols = [J[t] for t in range(k_) if t != c]
                        sub = np.array([[1 if rr == cc else 0 for cc in cols] for rr in rows], dtype=object).reshape(len(rows), len(cols))
                        cof = ((-1) ** (r + c)) * (minor(sub, list(range(len(rows))), list(range(len(cols)))) if rows else 1)
                        val += cof
            D[x, y] = val
    return D

DK = {k: {(a, b): Dk(k, a, b) for a in range(4) for b in range(4)} for k in range(1, 5)}

def rho(site, a, b, one=1):
    m = np.zeros((4, 4), dtype=object); m[:] = 0
    if site == 'A': m[a, b] = one
    else: m[b, a] = -one
    return m

def thetatilde(site, k, th, h):
    return th if site == 'A' else th + (k - 1) * h

def gfac(site, k, u, th, h):
    """scalar prefactor of the fused Lax (k=1:1; k=2: A u-th-h, F u-th; k=3: A (u-th-h)(u-th-2h), F (u-th)(u-th-h))"""
    if k <= 1: return 1
    if site == 'A':
        r = 1
        for j in range(1, k): r = r * (u - th - j * h)
        return r
    r = 1
    for j in range(0, k - 1): r = r * (u - th - j * h)
    return r

def gfac_reg(site, k, u, th, h, reg):
    if not reg: return gfac(site, k, u, th, h)
    js = range(1, k) if site == 'A' else range(0, k - 1)
    r = 1
    for j in js:
        f = u - th - j * h
        if f != 0: r = r * f
    return r

def lax_k(site, k, u, th, h):
    """fused Lax on wedge^k aux (x) site, as array [I,J,p,q] (site matrix acting on bras from the right), WITHOUT g prefactor"""
    S = subsets(k); n = len(S)
    L = np.zeros((n, n, 4, 4), dtype=object); L[:] = 0
    tt = thetatilde(site, k, th, h)
    for x in range(n):
        for p in range(4): L[x, x, p, p] = u - tt
    for a in range(4):
        for b in range(4):
            D = DK[k][(a, b)]; R = rho(site, a, b)
            for x in range(n):
                for y in range(n):
                    if D[x, y] != 0:
                        L[x, y] += h * D[x, y] * R
    return L

class Chain:
    def __init__(self, sites, th, h, G):
        self.sites = list(sites); self.th = list(th); self.h = h; self.G = G; self.N = len(sites)
        self.wG = {k: wedge(G, k) for k in range(1, 5)}
    def apply_t(self, k, u, psi, reg=False):
        """psi . t_k(u), psi shaped (4,)*N (object or numeric).  k=0 -> identity, k=4 -> qdet (scalar*psi).
        reg=True: drop the exactly-vanishing scalar prefactors (u-theta_s-jh)=0 ('hat' transfer matrix)."""
        if k == 0: return psi.copy()
        N = self.N; h = self.h
        n = len(subsets(k))
        pref = 1
        # Phi[I,J, sites...]
        Phi = np.empty((n, n) + psi.shape, dtype=psi.dtype)
        Phi[...] = 0
        for I in range(n): Phi[I, I] = psi
        for s in range(N):
            L = lax_k(self.sites[s], k, u, self.th[s], h)
            pref = pref * gfac_reg(self.sites[s], k, u, self.th[s], h, reg)
            if psi.dtype != object: L = L.astype(psi.dtype)
            # new Phi[I,J,...,q,...] = sum_K sum_p Phi[I,K,...,p,...] L[K,J,p,q]
            ax = 2 + s
            Phi = np.tensordot(Phi, L, axes=([1, ax], [0, 2]))  # -> [I, rest..., J, q]
            # move J to position 1 and q to position ax
            Phi = np.moveaxis(Phi, -2, 1)
            Phi = np.moveaxis(Phi, -1, ax)
        wg = self.wG[k]
        if psi.dtype != object: wg = wg.astype(psi.dtype)
        out = np.tensordot(wg, Phi, axes=([1, 0], [0, 1]))  # sum_{I,J} wg[J,I] Phi[I,J]
        # tensordot(wg,Phi,axes=([1,0],[0,1])) sums wg[a,b]*Phi[b,a]: a=J,b=I -> wg[J,I] Phi[I,J]  OK
        if psi.dtype != object and isinstance(pref, Fr): pref = float(pref)
        return pref * out

def prod_state(vecs, dtype=object):
    out = np.array(vecs[0], dtype=dtype)
    for v in vecs[1:]:
        out = np.multiply.outer(out, np.array(v, dtype=dtype))
    return out

def companion(z):
    from itertools import combinations
    def e(k):
        s = 0
        for c in combinations(z, k):
            p = 1
            for x in c: p = p * x
            s = s + p
        return s
    c = [e(k) for k in range(1, 5)]
    G = np.zeros((4, 4), dtype=object); G[:] = 0
    G[0, :] = [c[0], -c[1], c[2], -c[3]]
    G[1, 0] = 1; G[2, 1] = 1; G[3, 2] = 1
    return G

def diagm(z):
    G = np.zeros((4, 4), dtype=object); G[:] = 0
    for i in range(4): G[i, i] = z[i]
    return G

def matmul(A, B): return np.dot(A, B)

def inv(M):
    import sympy
    return np.array(sympy.Matrix(M.tolist()).inv().tolist(), dtype=object)

def det(M):
    import sympy
    return sympy.Matrix(M.tolist()).det()

def toFr(x):
    import sympy
    if isinstance(x, Fr): return x
    if isinstance(x, int): return Fr(x)
    x = sympy.nsimplify(x)
    return Fr(int(x.p), int(x.q))

def frM(M): return np.vectorize(toFr, otypes=[object])(np.array(M, dtype=object))

def site_op(psi, s, M):
    """psi . (M on site s)"""
    out = np.tensordot(psi, M, axes=([s], [0]))
    return np.moveaxis(out, -1, s)
