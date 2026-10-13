import sys, time, pickle, random, itertools
from fractions import Fraction
import numpy as np
from modlib import *
from modrun import mod_targets, localised_dim_mod, Wex_pair, vec_mod
from words import factor_list

def modmat_mul(A, B, p):
    """exact (A@B) mod p for float64 integer matrices, inner dim blocked by 128"""
    n = A.shape[1]
    out = np.zeros((A.shape[0], B.shape[1]))
    for b0 in range(0, n, 128):
        b1 = min(n, b0 + 128)
        out = np.mod(out + np.mod(A[:, b0:b1] @ B[b0:b1], p), p)
    return out

def mod_inv_matrix(A, p):
    n = A.shape[0]
    M = np.concatenate([np.mod(A, p), np.eye(n)], axis=1)
    for c in range(n):
        piv = None
        nz = np.nonzero(M[c:, c])[0]
        if len(nz) == 0: raise ZeroDivisionError('singular mod p')
        piv = c + nz[0]
        if piv != c: M[[c, piv]] = M[[piv, c]]
        M[c] = np.mod(M[c] * pow(int(M[c, c]), -1, p), p)
        col = M[:, c].copy(); col[c] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = np.mod(M[nzr] - np.mod(np.outer(col[nzr], M[c]), p), p)
    return M[:, n:]

class Cyclic:
    def __init__(self, L, thp, pi, seed=1, a0=Fraction(1, 3), verbose=True):
        self.L, self.thp = L, thp
        self.p = p = PRIMES[pi]
        th = {1: TH6, 2: TH6b, 3: TH6c}[thp]
        self.mc = mc = ModChain(L, th, H, Z4, p)
        v, u = Wex_pair(a0)
        self.psi0 = vec_mod([v, u] * (L // 2), p)
        rnd = random.Random(seed)
        u1, u2, u3 = [Fraction(rnd.randint(-97, 97), rnd.randint(11, 53)) for _ in range(3)]
        c2, c3 = rnd.randint(2, p - 2), rnd.randint(2, p - 2)
        def M(x):
            return np.mod(mc.apply(1, u1, x) + np.mod(c2 * mc.apply(2, u2, x), p) + np.mod(c3 * mc.apply(3, u3, x), p), p)
        basis = ModBasis(mc.dim, p)
        x = self.psi0.copy(); t0 = time.time()
        while True:
            if not basis.add(x[None, :])[0]: break
            x = M(x)
        self.basis = basis; self.r = basis.rank
        self.piv = np.array(basis.pivots)
        if verbose: print('cyclic dim', self.r, 'time', time.time() - t0, flush=True)
        self.F = factor_list(mc)
        self.T, self.w0 = mod_targets(L, a0, Z4, p)
        self.Tc = {k: vv[self.piv] for k, vv in self.T.items()}
        self.e0 = self.psi0[self.piv]
        self.ops = {}
    def coords(self, X):
        return X[..., self.piv]
    def op(self, f):
        """r x r matrix of factor f (dict) acting on coordinates (row-vector convention)."""
        key = (f['a'], f['s'], f['m'])
        if key not in self.ops:
            Y = np.concatenate([self.mc.apply(f['a'], f['u'], self.basis.rows[i0:i0 + 64].T).T for i0 in range(0, self.r, 64)], axis=0)
            Mm = Y[:, self.piv]
            # invariance check
            chk = modmat_mul(Mm, self.basis.rows, self.p)
            assert np.array_equal(chk, np.mod(Y, self.p)), 'not invariant'
            self.ops[key] = Mm
        return self.ops[key]
    def vec(self, coords):
        return modmat_mul(coords[None, :], self.basis.rows, self.p)[0]
    def loc_basis_coords(self):
        """coordinates of basis of S12 = span ∩ (C^16 x w0), mod p."""
        p = self.p; dimv = self.mc.dim; R = dimv // 16
        w0o = self.w0.astype(np.int64).astype(object)
        ww = int(np.dot(w0o, w0o) % p); winv = pow(ww, -1, p)
        P = ((np.eye(R, dtype=object) - np.outer(w0o, w0o) * winv) % p).astype(np.float64)
        rows = self.basis.rows.reshape(-1, 16, R)
        out = np.zeros_like(rows)
        for b0 in range(0, R, 128):
            b1 = min(R, b0 + 128)
            out = np.mod(out + np.mod(rows[:, :, b0:b1] @ P[b0:b1, :], p), p)
        flat = out.reshape(self.r, -1)     # r x (16R): map coords -> nonloc part: c |-> c @ flat
        return nullspace_left(flat, p)
    
def nullspace_left(A, p):
    """left nullspace {c : c @ A = 0 mod p} of r x m matrix A. returns k x r."""
    r, m = A.shape
    M = np.concatenate([np.mod(A, p), np.eye(r)], axis=1)
    row = 0
    for c in range(m):
        if row >= r: break
        nz = np.nonzero(M[row:, c])[0]
        if len(nz) == 0: continue
        piv = row + nz[0]
        if piv != row: M[[row, piv]] = M[[piv, row]]
        M[row] = np.mod(M[row] * pow(int(M[row, c]), -1, p), p)
        col = M[:, c].copy(); col[row] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = np.mod(M[nzr] - np.mod(np.outer(col[nzr], M[row]), p), p)
        row += 1
    return M[row:, m:]

def rank_mod(A, p):
    B = ModBasis(A.shape[1], p)
    fl = B.add(A.copy())
    return B.rank

def in_span(A, t, p):
    B = ModBasis(A.shape[1], p)
    B.add(A.copy())
    return B.contains(t)
