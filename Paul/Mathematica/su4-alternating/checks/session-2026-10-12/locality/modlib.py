"""Modular (mod p, p < 2^22, float64-exact BLAS) version of the transfer-matrix passes and rank/membership tests.
All rational inputs (theta, h, z, a) are reduced mod p; ranks are then exact up to accidental rank drops (prob ~ r/p).
"""
import numpy as np
from fractions import Fraction
from lib import *
import sympy

PRIMES = [int(sympy.prevprime(2 ** 22)), int(sympy.prevprime(2 ** 22 - 5000))]


def fm(x, p):
    x = Fraction(x)
    return (x.numerator % p) * pow(x.denominator % p, -1, p) % p


class ModChain:
    def __init__(self, nsites, th, h, z, p):
        self.p = p
        self.n = nsites
        self.types = ['A' if i % 2 == 0 else 'F' for i in range(nsites)]
        self.thf = th[:nsites]
        self.th = self.thf
        self.h = h
        self.hf = h
        self.zf = z
        self.dim = 4 ** nsites
        # exact-rational chain used only for building Lax tensors with object dtype, then reduce mod p
        G = np.zeros((4, 4), dtype=object)
        for i in range(4):
            for j in range(4):
                G[i, j] = Fraction(z[i]) if i == j else Fraction(0)
        self.G = G
        self.tw = {w: twist_matrix(w, G) for w in (1, 2, 3)}
        self._lax = {}

    def red(self, A):
        p = self.p
        out = np.zeros(A.shape, dtype=np.float64)
        it = np.nditer(np.empty(A.shape), flags=['multi_index'])
        for _ in it:
            out[it.multi_index] = fm(A[it.multi_index], p)
        return out

    def gprod(self, which, u):
        g = Fraction(1)
        h = Fraction(self.hf)
        for s in range(self.n):
            t = Fraction(self.thf[s])
            if which == 1:
                gs = 1
            elif which == 2:
                gs = (u - t - h) if self.types[s] == 'A' else (u - t)
            else:
                gs = (u - t - h) * (u - t - 2 * h) if self.types[s] == 'A' else (u - t) * (u - t - h)
            g = g * gs
        return g

    def laxes(self, which, u):
        key = (which, u)
        if key not in self._lax:
            if len(self._lax) > 300:
                self._lax.clear()
            Ls = []
            for s in range(self.n):
                Lx, g = lax_tensor(which, self.types[s], Fraction(u), Fraction(self.thf[s]), Fraction(self.hf), object)
                Ls.append(self.red(Lx))
            self._lax[key] = (Ls, fm(self.gprod(which, Fraction(u)), self.p))
        return self._lax[key]

    def apply(self, which, u, Psi, prefactor=True):
        """Psi float64 array (dim,B) of integers mod p."""
        p = self.p
        one_d = (Psi.ndim == 1)
        if one_d:
            Psi = Psi[:, None]
        Ls, gtot = self.laxes(which, u)
        d = 4 if which in (1, 3) else 6
        n = self.n
        B = Psi.shape[1]
        Phi = np.zeros((d, d, self.dim, B))
        for i in range(d):
            Phi[i, i] = Psi
        for s in range(n):
            Lx = Ls[s]
            Phi5 = Phi.reshape(d, d, 4 ** s, 4, 4 ** (n - s - 1) * B)
            new = np.tensordot(Phi5, Lx, axes=([1, 3], [0, 2]))
            new = np.mod(new, p)
            Phi = np.ascontiguousarray(np.transpose(new, (0, 3, 1, 4, 2))).reshape(d, d, self.dim, B)
        Gt = self.tw[which]
        out = np.zeros((self.dim, B))
        for i in range(d):
            for j in range(d):
                c = Gt[j, i]
                if c != 0:
                    out = np.mod(out + fm(c, p) * Phi[i, j], p)
        if prefactor:
            out = np.mod(out * gtot, p)
        return out[:, 0] if one_d else out


def mod_inv(a, p):
    return pow(int(a), -1, p)


class ModBasis:
    """incremental reduced basis (rows) mod p, float64 exact."""
    def __init__(self, dim, p):
        self.p = p
        self.dim = dim
        self.rows = np.zeros((0, dim))
        self.pivots = []   # pivot column per row (row normalised to 1 at pivot, zero at other pivots? not fully reduced)

    def _reduce(self, X):
        """reduce rows X (k x dim) against basis: X - X[:,pivots] @ rows (rows have unit pivot and are zero at others' pivots)."""
        p = self.p
        if len(self.pivots) == 0:
            return X
        piv = np.array(self.pivots)
        C = X[:, piv]
        # exact matmul in blocks of 128 inner terms
        out = X.copy()
        r = len(self.pivots)
        for b0 in range(0, r, 128):
            b1 = min(r, b0 + 128)
            out = np.mod(out - np.mod(C[:, b0:b1] @ self.rows[b0:b1], p), p)
        return out

    def add(self, X):
        """add rows X (k x dim); returns list of booleans: independent (was added) or not, in order."""
        p = self.p
        X = np.mod(X, p)
        flags = []
        X = self._reduce(X)
        # now reduce the chunk against itself sequentially
        for i in range(X.shape[0]):
            row = X[i]
            nzs = np.nonzero(row)[0]
            if len(nzs) == 0:
                flags.append(False)
                continue
            c = nzs[0]
            row = np.mod(row * mod_inv(row[c], p), p)
            # eliminate this pivot from existing basis rows and the remaining chunk rows
            if len(self.pivots):
                col = self.rows[:, c].copy()
                nzr = np.nonzero(col)[0]
                if len(nzr):
                    self.rows[nzr] = np.mod(self.rows[nzr] - np.mod(np.outer(col[nzr], row), p), p)
            if i + 1 < X.shape[0]:
                col = X[i + 1:, c].copy()
                nzr = np.nonzero(col)[0]
                if len(nzr):
                    X[i + 1 + nzr] = np.mod(X[i + 1 + nzr] - np.mod(np.outer(col[nzr], row), p), p)
            self.rows = np.vstack([self.rows, row[None, :]])
            self.pivots.append(c)
            flags.append(True)
        return flags

    def contains(self, x):
        x = np.mod(np.asarray(x, dtype=np.float64)[None, :], self.p)
        r = self._reduce(x)
        return not np.any(r)

    @property
    def rank(self):
        return len(self.pivots)

    def intersect_dim_localised(self, w0_mod):
        """dimension of basis-span intersect (C^16 x w0): solve; uses nullspace of the map c -> proj_{w0^perp}(psi)"""
        raise NotImplementedError
