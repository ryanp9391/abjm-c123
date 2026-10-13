"""Taylor coefficients of t_a(u0 + eta) in eta (exact mod p) and series words."""
import numpy as np
from fractions import Fraction as Fr
from chain import *
from words import *


def frac_inv(V):
    n = len(V)
    A = [list(map(Fr, row)) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(V)]
    for c in range(n):
        pr = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[pr] = A[pr], A[c]
        iv = 1 / A[c][c]
        A[c] = [x * iv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def lagrange_taylor_weights(nodes, order, p=P):
    """weights w[m][i] with sum_i w[m][i] f(nodes[i]) = f^{(m)}(0)/m! for polys of degree < len(nodes)."""
    n = len(nodes)
    V = [[Fr(x) ** k for k in range(n)] for x in nodes]
    Vi = frac_inv(V)
    F = ModP(p)
    return [[F.c(Vi[m][i]) if m < n else 0 for i in range(n)] for m in range(order + 1)]


def t_taylor(ch, a, u0, order):
    """list of (order+1) matrices: Taylor coeffs of t_a(u0+eta)."""
    F = ch.F
    if a == 0:
        return [F.eye(ch.N)] + [F.zeros(ch.N) for _ in range(order)]
    deg = ch.n * a
    nodes = [Fr(s) for s in range(-(deg // 2), deg - deg // 2 + 1)]
    nodes = nodes[:deg + 1]
    W = lagrange_taylor_weights(nodes, order, F.p)
    vals = []
    for x in nodes:
        vals.append(ch.t(a, Fr(u0) + x))
        ch._tcache.pop((a, Fr(u0) + x), None)
        ch._Tcache.clear()
    out = []
    for m in range(order + 1):
        acc = F.zeros(ch.N)
        for i, x in enumerate(nodes):
            acc = (acc + W[m][i] * vals[i]) % F.p
        out.append(acc)
    return out


# ---- truncated series of matrices / vectors: arrays shape (N+1, ...) ----
def ser_mm(A, B, p=P):
    N = A.shape[0] - 1
    out = np.zeros((N + 1,) + A.shape[1:-1] + B.shape[2:], dtype=np.int64)
    for i in range(N + 1):
        for j in range(N + 1 - i):
            out[i + j] = (out[i + j] + (A[i] @ B[j]) % p) % p
    return out


def ser_vm(v, B, p=P):
    """v: (N+1, n) ; B: (N+1, n, n)"""
    N = v.shape[0] - 1
    out = np.zeros_like(v)
    for i in range(N + 1):
        for j in range(N + 1 - i):
            out[i + j] = (out[i + j] + (v[i] @ B[j]) % p) % p
    return out


def ser_scale_shift(coeffs, lam, p=P):
    """given Taylor coeffs C_m of f(u0+eta), return series of f(u0 + lam*eps) in eps."""
    lam = ModP(p).c(lam)
    return np.array([(pow(lam, m, p) * C) % p for m, C in enumerate(coeffs)])


def ser_words(ch, ops):
    """ops: dict with 'A1','A3' : list over a of series (N+1,n,n); 'T2','T4': series.
    returns (N+1, 256, n) array of words."""
    p = ch.F.p
    N = ops['T2'].shape[0] - 1
    vac = vacuum(ch.types, ch.F)
    v0 = np.zeros((N + 1, ch.N), dtype=np.int64)
    v0[0] = vac
    rows = {}
    for a1 in range(4):
        v = ser_vm(v0, ops['A1'][a1], p)
        for k2 in range(4):
            if k2 > 0:
                v = ser_vm(v, ops['T2'], p)
            for a3 in range(4):
                w = ser_vm(v, ops['A3'][a3], p)
                for k4 in range(4):
                    if k4 > 0:
                        w = ser_vm(w, ops['T4'], p)
                    rows[(a1, k2, a3, k4)] = w
    return np.stack([rows[k] for k in WIDX], axis=1)  # (N+1, 256, n)


def valuation_basis(Ms, p=P):
    """Ms: (N+1, r, n) series rows w_i(eps). Gaussian elimination over GF(p)[[eps]] (truncated):
    returns list of (n_j, combination series S_j (N+1, r), leading vector v_j (n,)) such that
    sum_i S_j[i](eps) w_i(eps) = eps^{n_j} (v_j + O(eps)), v_j independent.
    Strategy: process orders m = 0..N: at order m collect leading vectors of remaining rows that are
    independent of previously selected leading vectors; rows whose leading part is dependent get reduced
    (subtract combination) and their valuation increases."""
    N = Ms.shape[0] - 1
    r, n = Ms.shape[1], Ms.shape[2]
    # each row: series vector (N+1, n) and combination (N+1, r)
    rows = [Ms[:, i, :].copy() for i in range(r)]
    combs = []
    for i in range(r):
        c = np.zeros((N + 1, r), dtype=np.int64)
        c[0, i] = 1
        combs.append(c)
    vals = []
    for i in range(r):
        nzm = [m for m in range(N + 1) if np.any(rows[i][m])]
        vals.append(nzm[0] if nzm else N + 1)
    selected = []  # (val, comb, lead)
    # pivot structure for selected leading vectors: rref kept incrementally
    piv_rows = []  # list of (pivot col, normalized vector, index into selected)
    active = list(range(r))
    for m in range(N + 1):
        changed = True
        while changed:
            changed = False
            for i in list(active):
                if vals[i] != m:
                    continue
                lead = rows[i][m].copy()
                # reduce lead against pivots, record combination
                coef = []
                for (pc, vec, sidx) in piv_rows:
                    c = lead[pc]
                    if c:
                        lead = (lead - c * vec) % p
                        coef.append((c, sidx))
                if np.any(lead):
                    # new independent leading vector
                    nzc = np.nonzero(lead)[0][0]
                    inv = pow(int(lead[nzc]), p - 2, p)
                    vec = (lead * inv) % p
                    # keep rref-ish: reduce existing pivot vectors by new one not needed for correctness
                    selected.append([m, combs[i].copy(), rows[i][m].copy(), rows[i].copy()])
                    # store reduced vec referencing selected index; to subtract we need the full series of the
                    # selected row times eps^(m - m_sel) ... handled below via coefficients on selected rows
                    piv_rows.append((nzc, vec, len(selected) - 1, ) if False else (nzc, vec, len(selected) - 1))
                    # express vec in terms of selected leading vectors: store transformation
                    active.remove(i)
                else:
                    # dependent: subtract sum_c c * eps^{m - m_s} * (row_s / leadnorm) ... need exact representation
                    # recompute: find coefficients x_s with lead_i = sum_s x_s lead_s (lead_s = leading vec of selected s)
                    L = np.array([s[2] for s in selected], dtype=np.int64)
                    x = solve_in_span(L, rows[i][m], p)
                    newrow = rows[i].copy()
                    newcomb = combs[i].copy()
                    for sidx, xs in enumerate(x):
                        if xs == 0:
                            continue
                        ms = selected[sidx][0]
                        sh = m - ms
                        srow = selected[sidx][3]
                        scomb = selected[sidx][1]
                        # subtract xs * eps^sh * srow
                        newrow[sh:] = (newrow[sh:] - xs * srow[:N + 1 - sh]) % p
                        newcomb[sh:] = (newcomb[sh:] - xs * scomb[:N + 1 - sh]) % p
                    rows[i] = newrow
                    combs[i] = newcomb
                    nzm = [mm for mm in range(N + 1) if np.any(newrow[mm])]
                    vals[i] = nzm[0] if nzm else N + 1
                    changed = True
    return selected, [(i, vals[i]) for i in active]


def solve_in_span(L, v, p=P):
    """x with x . L = v (L rows independent); raises if not in span."""
    k = L.shape[0]
    A = np.concatenate([L.T % p, (v % p).reshape(-1, 1)], axis=1)
    R, piv = rref_mod(A, p)
    assert k not in piv, 'not in span'
    x = np.zeros(k, dtype=np.int64)
    for row, pc in enumerate(piv):
        x[pc] = R[row, k]
    return x
