# fast bra-action of twisted transfer matrices for an L-site A/F chain; ring = mod p (exact) or float
import numpy as np
from fractions import Fraction
from itertools import permutations, combinations
P = 2147483647
def perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s
class Ring:
    def __init__(self, mode): self.mode = mode
    def c(self, x):  # convert Fraction/int to scalar
        x = Fraction(x)
        if self.mode == 'mod': return (x.numerator % P) * pow(x.denominator % P, P - 2, P) % P
        return float(x)
    def red(self, a): return a % P if self.mode == 'mod' else a
    def zeros(self, n): return np.zeros(n, dtype=np.int64 if self.mode == 'mod' else np.float64)
class FastChain:
    def __init__(self, types, th, h, cvec, ring):
        self.types = types; self.L = len(types); self.R = ring
        self.th = [ring.c(t) for t in th]; self.h = ring.c(h)
        c1, c2, c3, c4 = [ring.c(x) for x in cvec]
        neg = (lambda x: (-x) % P) if ring.mode == 'mod' else (lambda x: -x)
        self.G = [[c1, neg(c2), c3, neg(c4)], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
        if ring.mode == 'float': self.G = [[float(x) for x in row] for row in self.G]
        self.dim = 4**self.L
    def local(self, psi, s, a, b, u):
        # psi . L_s(u)_{ab}, psi flat array
        R = self.R
        t = psi.reshape((4,) * self.L)
        out = R.zeros(t.shape)
        if a == b: out = R.red(t * R.red(u - self.th[s]))
        # h rho_s(E_ab)
        idx_src = [slice(None)] * self.L; idx_dst = [slice(None)] * self.L
        if self.types[s] == 'A':   # psi e_ab : component b <- psi component a
            idx_src[s] = a; idx_dst[s] = b; sgn = 1
        else:                      # psi (-e_ba): component a <- -psi component b
            idx_src[s] = b; idx_dst[s] = a; sgn = -1
        add = R.red(t[tuple(idx_src)] * (self.h if sgn == 1 else R.red(-self.h)))
        out[tuple(idx_dst)] = R.red(out[tuple(idx_dst)] + add)
        return out.reshape(-1)
    def braTG(self, psi, u):
        R = self.R
        # X[i][k] = psi . L1_{ik} ... Ls_{.k}
        X = [[self.local(psi, 0, i, k, u) for k in range(4)] for i in range(4)]
        for s in range(1, self.L):
            Y = [[R.zeros(self.dim) for k in range(4)] for i in range(4)]
            for i in range(4):
                for m in range(4):
                    for k in range(4):
                        Y[i][k] = R.red(Y[i][k] + self.local(X[i][m], s, m, k, u))
            X = Y
        TG = [[R.red(sum(R.red(X[i][l] * self.G[l][j]) for l in range(4))) for j in range(4)] for i in range(4)]
        return TG
    def bra_qmin(self, psi, I, J, u):
        R = self.R
        out = R.zeros(self.dim)
        n = len(J)
        for p in permutations(range(n)):
            sg = perm_sign(p)
            v = psi
            for k in range(n):
                uk = R.red(u - k * self.h)
                v = self.braTG(v, uk)[I[p[k]]][J[k]]
            out = R.red(out + (sg * v if R.mode == 'float' else (v if sg == 1 else (-v) % P)))
        return out
    def bra_t(self, psi, a, u):
        if a == 0: return psi
        R = self.R
        if a == 1:
            TG = self.braTG(psi, u)
            return R.red(sum(TG[i][i] for i in range(4)))
        out = R.zeros(self.dim)
        for I in combinations(range(4), a):
            out = R.red(out + self.bra_qmin(psi, I, I, u))
        return out
def rank_mod(M):
    M = np.array(M, dtype=object) % P
    M = M.copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c] % P != 0), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), P - 2, P)
        M[r] = (M[r] * inv) % P
        for i in range(rows):
            if i != r and M[i, c] != 0:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        r += 1
        if r == rows: break
    return r
def inv_mod(M):
    n = M.shape[0]
    A = np.concatenate([np.array(M, dtype=object) % P, np.eye(n, dtype=object)], axis=1)
    for c in range(n):
        piv = next(i for i in range(c, n) if A[i, c] % P != 0)
        A[[c, piv]] = A[[piv, c]]
        inv = pow(int(A[c, c]), P - 2, P)
        A[c] = (A[c] * inv) % P
        for i in range(n):
            if i != c and A[i, c] != 0:
                A[i] = (A[i] - A[i, c] * A[c]) % P
    return A[:, n:]
