# one-pair (A F) and two-pair chains: Lax, twisted monodromy, transfer matrices acting on bras (row vectors)
import sympy as sp
from itertools import permutations, combinations
from sympy.combinatorics import Permutation

def E(i, j, n=4):
    m = sp.zeros(n); m[i, j] = 1; return m

def siteop(L, s, M):
    # operator M on site s (0-based) of an L-site chain, as dense sympy matrix
    mats = [M if b == s else sp.eye(4) for b in range(L)]
    out = mats[0]
    for m in mats[1:]:
        out = sp.kronecker_product(out, m)
    return out

def perm_sign(p):
    return Permutation(list(p)).signature()

class Chain:
    def __init__(self, types, th, h, G):
        self.types = types; self.L = len(types); self.th = th; self.h = h; self.G = G
        self.dim = 4**self.L
        # site generators rho_s(E_ij) as dense matrices
        self.rho = {}
        for s, t in enumerate(types):
            for i in range(4):
                for j in range(4):
                    M = E(i, j) if t == 'A' else -E(j, i)
                    self.rho[(s, i, j)] = siteop(self.L, s, M)
        self.cache = {}
    def lax(self, s, u):
        return [[(u - self.th[s]) * (sp.eye(self.dim) if i == j else sp.zeros(self.dim)) + self.h * self.rho[(s, i, j)]
                 for j in range(4)] for i in range(4)]
    def TG(self, u):
        key = sp.simplify(u)
        if key in self.cache: return self.cache[key]
        T = self.lax(0, u)
        for s in range(1, self.L):
            Ls = self.lax(s, u)
            T = [[sum((T[i][k] * Ls[k][j] for k in range(4)), sp.zeros(self.dim)) for j in range(4)] for i in range(4)]
        TGm = [[sum((T[i][l] * self.G[l, j] for l in range(4)), sp.zeros(self.dim)) for j in range(4)] for i in range(4)]
        TGm = [[TGm[i][j].applyfunc(sp.expand) for j in range(4)] for i in range(4)]
        self.cache[key] = TGm
        return TGm
    def bra_qmin(self, psi, I, J, u):
        # <psi| qmin_{I,J}(u) ; rows antisymmetrised, column k at u-(k-1)h
        n = len(J)
        mats = [self.TG(u - k * self.h) for k in range(n)]
        out = sp.zeros(1, self.dim)
        for p in permutations(range(n)):
            sg = perm_sign(p)
            v = psi
            for k in range(n):
                v = v * mats[k][I[p[k]]][J[k]]
            out += sg * v
        return out.applyfunc(sp.expand)
    def bra_t(self, psi, a, u):
        if a == 0: return psi
        out = sp.zeros(1, self.dim)
        for I in combinations(range(4), a):
            out += self.bra_qmin(psi, I, I, u)
        return out.applyfunc(sp.expand)
    def op_t(self, a, u):
        # full operator matrix t_a(u) (rows = bra action): build by acting on unit bras
        rows = []
        for n in range(self.dim):
            e = sp.zeros(1, self.dim); e[0, n] = 1
            rows.append(self.bra_t(e, a, u))
        return sp.Matrix.vstack(*rows)

def companion(c):
    c1, c2, c3, c4 = c
    return sp.Matrix([[c1, -c2, c3, -c4], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])

def esym(z):
    return [sum(sp.prod(s) for s in combinations(z, k)) for k in range(1, 5)]
