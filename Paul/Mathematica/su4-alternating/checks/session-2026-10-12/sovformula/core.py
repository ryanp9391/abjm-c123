"""Core: 4-site A F A F su(4) chain, companion/diagonal twist, t1,t2,t3, SoV words.
Two arithmetic modes: 'float' (numpy float64/complex128) and 'mod' (exact mod prime p, int64)."""
import numpy as np
from fractions import Fraction as Fr
from itertools import permutations, combinations, product
import math

P1 = 67108859
P2 = 67108837


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


class Ring:
    """float or mod-p arithmetic on numpy arrays"""

    def __init__(self, mode='float', p=P1):
        self.mode = mode
        self.p = p
        self.dtype = np.int64 if mode == 'mod' else np.float64

    def c(self, x):  # scalar Fraction/int -> ring element
        if self.mode == 'mod':
            x = Fr(x)
            return (x.numerator % self.p) * pow(x.denominator % self.p, -1, self.p) % self.p
        if isinstance(x, complex):
            return x
        return float(x)

    def arr(self, a):
        a = np.array(a, dtype=object)
        if self.mode == 'mod':
            return np.vectorize(self.c, otypes=[np.int64])(a) if a.size else a.astype(np.int64)
        out = np.vectorize(lambda x: complex(x) if isinstance(x, complex) else float(x), otypes=[complex])(a)
        if np.all(out.imag == 0):
            return out.real.astype(float)
        return out

    def red(self, a):
        return a % self.p if self.mode == 'mod' else a

    def mm(self, a, b):
        if self.mode == 'mod':
            p = self.p
            lo = b & 0x1FFF
            hi = b >> 13
            r1 = (a @ lo) % p
            r2 = (a @ hi) % p
            return (r1 + (r2 * 8192) % p) % p
        return a @ b

    def kron(self, a, b):
        return self.red(np.kron(a, b))

    def smul(self, x, a):  # scalar times array
        if self.mode == 'mod':
            return (self.c(x) * a) % self.p
        return self.c(x) * a

    def inv(self, A):
        if self.mode == 'float':
            return np.linalg.inv(A)
        return modinv(A, self.p)

    def rank(self, A):
        if self.mode == 'float':
            return np.linalg.matrix_rank(A)
        return modrank(A, self.p)


def rref_mod(A, p):
    A = A.copy() % p
    n, m = A.shape
    r = 0
    piv = []
    for c in range(m):
        pv = None
        col = A[r:, c]
        nzs = np.nonzero(col)[0]
        if len(nzs) == 0:
            continue
        pv = r + nzs[0]
        if pv != r:
            A[[r, pv]] = A[[pv, r]]
        inv = pow(int(A[r, c]), -1, p)
        A[r] = (A[r] * inv) % p
        f = A[:, c].copy()
        f[r] = 0
        nz = np.nonzero(f)[0]
        if len(nz):
            # avoid overflow: f<p<2^26, A<2^26 -> product < 2^52 ok
            A[nz] = (A[nz] - (f[nz, None] * A[r][None, :]) % p) % p
        piv.append(c)
        r += 1
        if r == n:
            break
    return A[:r], piv


def modrank(A, p):
    return len(rref_mod(A, p)[1])


def modinv(A, p):
    n = A.shape[0]
    Aug = np.concatenate([A % p, np.eye(n, dtype=np.int64)], axis=1)
    Rr, piv = rref_mod(Aug, p)
    if piv[:n] != list(range(n)):
        raise ValueError("singular mod p")
    return Rr[:, n:]


def nullspace_mod(A, p):
    """right nullspace basis (columns) of A mod p"""
    Rr, piv = rref_mod(A, p)
    n, m = A.shape
    free = [c for c in range(m) if c not in piv]
    N = np.zeros((m, len(free)), dtype=np.int64)
    for k, fc in enumerate(free):
        N[fc, k] = 1
        for i, pc in enumerate(piv):
            N[pc, k] = (-Rr[i, fc]) % p
    return N


def ratrec(a, p):
    """rational reconstruction of a mod p"""
    a = int(a) % p
    r0, r1 = p, a
    s0, s1 = 0, 1
    B = int(math.isqrt(p // 2))
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


# ---------------- parameters ----------------
h0 = Fr(1, 2)
th0 = [Fr(2, 7), Fr(-3, 5), Fr(8, 11), Fr(-11, 8)]
z0 = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
types = ['A', 'F', 'A', 'F']


def esym(z):
    e1 = sum(z)
    e2 = sum(z[i] * z[j] for i, j in combinations(range(4), 2))
    e3 = sum(z[i] * z[j] * z[k] for i, j, k in combinations(range(4), 3))
    e4 = z[0] * z[1] * z[2] * z[3]
    return [e1, e2, e3, e4]


def Gc_of(z):
    c = esym(z)
    return [[c[0], -c[1], c[2], -c[3]], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]


def Vand(z):
    return [[z[j] ** 3 for j in range(4)], [z[j] ** 2 for j in range(4)], [z[j] for j in range(4)], [1] * 4]


class Chain:
    def __init__(self, R, th=th0, h=h0, G=None, z=z0):
        self.R = R
        self.th = th
        self.h = h
        self.z = z
        self.G = G if G is not None else Gc_of(z)
        self.Gr = R.arr(self.G)
        self._cache = {}

    def lax(self, s, u):
        R = self.R
        ty = types[s]
        d = u - self.th[s]
        L = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                m = [[0] * 4 for _ in range(4)]
                if i == j:
                    for k in range(4):
                        m[k][k] = d
                if ty == 'A':
                    m[i][j] = m[i][j] + self.h
                else:
                    m[j][i] = m[j][i] - self.h
                L[i][j] = R.arr(m)
        return L

    def mono(self, u):
        R = self.R
        L = [self.lax(s, u) for s in range(4)]

        def pairprod(La, Lb):
            P = [[None] * 4 for _ in range(4)]
            for i in range(4):
                for l in range(4):
                    acc = None
                    for k in range(4):
                        t = R.kron(La[i][k], Lb[k][l])
                        acc = t if acc is None else R.red(acc + t)
                    P[i][l] = acc
            return P
        P12 = pairprod(L[0], L[1])
        P34 = pairprod(L[2], L[3])
        T = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                acc = None
                for l in range(4):
                    t = R.kron(P12[i][l], P34[l][j])
                    acc = t if acc is None else R.red(acc + t)
                T[i][j] = acc
        return T

    def monoG(self, u):
        key = ('TG', u)
        if key in self._cache:
            return self._cache[key]
        R = self.R
        T = self.mono(u)
        Gr = self.Gr
        TG = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                acc = np.zeros_like(T[0][0])
                if R.mode == 'float' and np.iscomplexobj(Gr):
                    acc = acc.astype(complex)
                for k in range(4):
                    g = Gr[k, j]
                    if g != 0:
                        acc = R.red(acc + T[i][k] * g)
                TG[i][j] = acc
        self._cache[key] = TG
        return TG

    def t(self, a, u):
        key = ('t', a, u)
        if key in self._cache:
            return self._cache[key]
        R = self.R
        if a == 0:
            res = R.arr(np.eye(256, dtype=int))
        else:
            TGs = [self.monoG(u - k * self.h) for k in range(a)]
            res = None
            for I in combinations(range(4), a):
                for sg in permutations(range(a)):
                    sgn = perm_sign(sg)
                    prod = TGs[0][I[sg[0]]][I[0]]
                    for k in range(1, a):
                        prod = R.mm(prod, TGs[k][I[sg[k]]][I[k]])
                    term = prod if sgn > 0 else R.red(-prod)
                    res = term if res is None else R.red(res + term)
        self._cache[key] = res
        return res


def kronvecs(R, vs):
    b = R.arr(vs[0]) if not isinstance(vs[0], np.ndarray) else vs[0]
    for v in vs[1:]:
        vv = R.arr(v) if not isinstance(v, np.ndarray) else v
        b = R.red(np.kron(b, vv))
    return b


def e(i):
    v = [0] * 4
    v[i - 1] = 1
    return v


wordsIdx = list(product(range(4), repeat=4))


def sov_words(ch):
    R = ch.R
    th = ch.th
    h = ch.h
    vac = kronvecs(R, [e(1), e(4), e(1), e(4)])
    A1 = [R.mm(vac[None, :], ch.t(a, th[0] - h))[0] for a in range(4)]
    T2 = ch.t(1, th[1])
    T4 = ch.t(1, th[3])
    T3 = [ch.t(a, th[2] - h) for a in range(4)]
    rows = []
    for (a1, k2, a3, k4) in wordsIdx:
        v = A1[a1][None, :]
        for _ in range(k2):
            v = R.mm(v, T2)
        v = R.mm(v, T3[a3])
        for _ in range(k4):
            v = R.mm(v, T4)
        rows.append(v[0])
    return np.array(rows), vac
