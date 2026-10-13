"""exact linear algebra mod P (eng.P) for tall systems A c = b (A: D x m, small m)."""
import numpy as np
from eng import P

def solve(cols, b):
    """cols: list of m vectors (int64 mod P, any shape), b: vector. Returns (rank, consistent, c (with free vars 0), nullspace basis list)."""
    m = len(cols)
    A = np.stack([c.reshape(-1) % P for c in cols] + [b.reshape(-1) % P], axis=1).astype(np.int64)
    D = A.shape[0]
    piv_cols = []; r = 0
    for c in range(m + 1):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), P - 2, P)
        A[r] = (A[r] * inv) % P
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % P) % P
        piv_cols.append(c); r += 1
        if r == D: break
    consistent = m not in piv_cols
    rank = len([c for c in piv_cols if c < m])
    sol = np.zeros(m, dtype=np.int64)
    for i, c in enumerate(piv_cols):
        if c < m: sol[c] = A[i, m]
    free = [c for c in range(m) if c not in piv_cols]
    null = []
    for f in free:
        v = np.zeros(m, dtype=np.int64); v[f] = 1
        for i, c in enumerate(piv_cols):
            if c < m: v[c] = (-A[i, f]) % P
        null.append(v)
    return rank, consistent, sol, null

def ratrec(x, p=P):
    """rational reconstruction of x mod p -> Fraction (Wang), bound sqrt(p/2)"""
    from fractions import Fraction
    import math
    x %= p
    B = int(math.isqrt(p // 2))
    r0, r1 = p, x; s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)
