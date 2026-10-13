"""su(4) alternating chain A F A F ..., transfer matrices from quantum minors.
Fields: ModP (GF(p), int64 numpy), Flt (complex128 numpy).
Conventions: Lax L_s(u)_ij = (u-th_s) d_ij + h rho_s(E_ij), rho_A(E_ij)=e_ij, rho_F(E_ij)=-e_ji,
T = L_1 ... L_n (site 1 leftmost, most significant tensor index), T^G = T.G,
qmin_{I,J}(u) = sum_sigma sgn T^G_{I_s(1),J_1}(u) T^G_{I_s(2),J_2}(u-h)...,
t_a(u) = sum_{|I|=a} qmin_{I,I}(u). Bras are row vectors, psi . M."""
import numpy as np
from fractions import Fraction as Fr
from itertools import permutations, combinations

P = 16777213  # prime < 2^24


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


class ModP:
    def __init__(self, p=P):
        self.p = p

    def c(self, x):
        x = Fr(x)
        return (x.numerator % self.p) * pow(x.denominator % self.p, self.p - 2, self.p) % self.p

    def arr(self, a):
        return np.asarray(a, dtype=np.int64) % self.p

    def mm(self, A, B):
        return (A @ B) % self.p

    def red(self, A):
        return A % self.p

    def zeros(self, n):
        return np.zeros((n, n), dtype=np.int64)

    def eye(self, n):
        return np.eye(n, dtype=np.int64)

    def scal(self, s, A):
        return (int(s) * A) % self.p

    def kron(self, A, B):
        return np.kron(A, B) % self.p


class Flt:
    def __init__(self, dtype=np.complex128):
        self.dtype = dtype

    def c(self, x):
        return complex(x) if self.dtype == np.complex128 else float(x)

    def arr(self, a):
        return np.asarray(a, dtype=self.dtype)

    def mm(self, A, B):
        return A @ B

    def red(self, A):
        return A

    def zeros(self, n):
        return np.zeros((n, n), dtype=self.dtype)

    def eye(self, n):
        return np.eye(n, dtype=self.dtype)

    def scal(self, s, A):
        return s * A

    def kron(self, A, B):
        return np.kron(A, B)


def e(i, j, n=4):
    m = np.zeros((n, n), dtype=np.int64)
    m[i, j] = 1
    return m


def esym(z, k):
    s = Fr(0)
    for c in combinations(z, k):
        pr = Fr(1)
        for x in c:
            pr *= Fr(x)
        s += pr
    return s


def companion(z):
    c = [esym(z, k) for k in range(1, 5)]
    return [[c[0], -c[1], c[2], -c[3]], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]


class Chain:
    def __init__(self, field, types, th, h, G):
        """th, h: rationals; G: 4x4 rationals. field ModP or Flt"""
        self.F = field
        self.types = types
        self.n = len(types)
        self.N = 4 ** self.n
        self.th = [Fr(x) for x in th]
        self.h = Fr(h)
        self.G = [[Fr(x) for x in r] for r in G]
        self._Tcache = {}
        self._tcache = {}
        self.rho = {}
        for ty in ('A', 'F'):
            for i in range(4):
                for j in range(4):
                    self.rho[(ty, i, j)] = e(i, j) if ty == 'A' else -e(j, i)

    def lax_entry(self, s, u, i, j):
        F = self.F
        ty = self.types[s]
        if isinstance(F, ModP):
            m = (self.rho[(ty, i, j)] * F.c(self.h)) % F.p
            if i == j:
                m = (m + F.c(u - self.th[s]) * np.eye(4, dtype=np.int64)) % F.p
            return m
        m = F.arr(self.rho[(ty, i, j)]) * F.c(self.h)
        if i == j:
            m = m + F.c(u - self.th[s]) * np.eye(4)
        return m

    def T(self, u):
        u = Fr(u)
        if u in self._Tcache:
            return self._Tcache[u]
        F = self.F
        n = self.n
        Ls = [[[self.lax_entry(s, u, i, j) for j in range(4)] for i in range(4)] for s in range(n)]
        cur = [[Ls[0][i][j] for j in range(4)] for i in range(4)]
        for s in range(1, n):
            new = [[None] * 4 for _ in range(4)]
            for i in range(4):
                for j in range(4):
                    acc = None
                    for k in range(4):
                        term = F.kron(cur[i][k], Ls[s][k][j])
                        acc = term if acc is None else F.red(acc + term)
                    new[i][j] = acc
            cur = new
        TG = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                acc = F.zeros(self.N)
                for k in range(4):
                    g = self.G[k][j]
                    if g != 0:
                        acc = F.red(acc + F.scal(F.c(g), cur[i][k]))
                TG[i][j] = acc
        self._Tcache[u] = (cur, TG)
        return cur, TG

    def qmin(self, I, J, u, twisted=True):
        u = Fr(u)
        F = self.F
        a = len(I)
        Ts = [self.T(u - k * self.h)[1 if twisted else 0] for k in range(a)]
        acc = F.zeros(self.N)
        for sg in permutations(range(a)):
            sgn = perm_sign(sg)
            M = Ts[0][I[sg[0]]][J[0]]
            for k in range(1, a):
                M = F.mm(M, Ts[k][I[sg[k]]][J[k]])
            acc = F.red(acc + (M if sgn > 0 else -M))
        return acc

    def t(self, a, u):
        u = Fr(u)
        key = (a, u)
        if key in self._tcache:
            return self._tcache[key]
        F = self.F
        if a == 0:
            r = F.eye(self.N)
        else:
            r = F.zeros(self.N)
            for I in combinations(range(4), a):
                r = F.red(r + self.qmin(I, I, u))
        self._tcache[key] = r
        return r


def vacuum(types, F):
    v = np.array([1], dtype=np.int64)
    for ty in types:
        ev = np.zeros(4, dtype=np.int64)
        ev[0 if ty == 'A' else 3] = 1
        v = np.kron(v, ev)
    if isinstance(F, ModP):
        return v
    return F.arr(v)


def rank_mod(M, p=P):
    A = np.array(M, dtype=np.int64) % p
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        nzr = np.nonzero(A[r:, c])[0]
        if len(nzr) == 0:
            continue
        piv = r + nzr[0]
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nz = np.nonzero(col)[0]
        if len(nz):
            A[nz] = (A[nz] - np.outer(col[nz], A[r])) % p
        r += 1
        if r == rows:
            break
    return r


def rref_mod(M, p=P):
    A = np.array(M, dtype=np.int64) % p
    rows, cols = A.shape
    r = 0
    piv = []
    for c in range(cols):
        if r == rows:
            break
        nzr = np.nonzero(A[r:, c])[0]
        if len(nzr) == 0:
            continue
        pr = r + nzr[0]
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nz = np.nonzero(col)[0]
        if len(nz):
            A[nz] = (A[nz] - np.outer(col[nz], A[r])) % p
        piv.append(c)
        r += 1
    return A[:r], piv


def left_null_mod(M, p=P):
    """basis of {c : c . M = 0} (rows of result)."""
    At = np.array(M, dtype=np.int64).T % p  # cols x rows ; solve At c = 0
    R, piv = rref_mod(At, p)
    ncol = At.shape[1]
    free = [c for c in range(ncol) if c not in piv]
    out = []
    for f in free:
        v = np.zeros(ncol, dtype=np.int64)
        v[f] = 1
        for k, pc in enumerate(piv):
            v[pc] = (-R[k, f]) % p
        out.append(v)
    return np.array(out, dtype=np.int64).reshape(len(out), ncol)


def solve_left_mod(M, b, p=P):
    """c with c . M = b (M square invertible mod p)."""
    n = M.shape[0]
    A = np.concatenate([np.array(M, dtype=np.int64).T % p, (np.array(b, dtype=np.int64) % p).reshape(-1, 1)], axis=1)
    R, piv = rref_mod(A, p)
    assert piv[:n] == list(range(n)), 'singular'
    return R[:n, n]


def ratrec(a, p=P):
    a = int(a) % p
    r0, r1 = p, a
    s0, s1 = 0, 1
    bound = int((p / 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)
