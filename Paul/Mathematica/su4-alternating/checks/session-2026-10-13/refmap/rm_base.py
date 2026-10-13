"""rm_base.py -- pair-level machinery for the K_kappa reference-admissibility map (python FALLBACK, written because the shared
wolfbook kernel was monopolised by another agent for >10 min on 2026-10-13; the same computations are in
reference-admissibility-Kkappa.wb).  Exact Fractions.  Conventions = pgbase.py (A F, bras row vectors, psi -> psi.M).
One-pair chain: sites (A, F), theta = (2h, 0) = (1, 0), h = 1/2, diagonal twist Lambda = diag(zz)."""
import sys, os, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'session-2026-10-12', 'pairgroup', 'py'))
import numpy as np
from pgbase import Chain, diagm, prod_state

h = Fr(1, 2)

def F(x): return Fr(x)

def vec(*xs): return [Fr(x) for x in xs]

_OPS = {}
def pair_ops_z(zz):
    key = tuple(zz)
    if key not in _OPS: _OPS[key] = pair_ops(diagm(list(zz)))
    return _OPS[key]

def pair_ops(G, th=(Fr(1), Fr(0))):
    """16x16 matrices (rows = images of basis bras e_a (x) e_b) of X1 = t1(thF), X2 = t2(thF+h), X3 = t3hat(thF+h)."""
    ch = Chain(['A', 'F'], list(th), h, G)
    thF = th[1]
    specs = [(1, thF, False), (2, thF + h, False), (3, thF + h, True)]
    out = []
    for k, u, reg in specs:
        M = [[Fr(0)] * 16 for _ in range(16)]
        for a in range(4):
            for b in range(4):
                psi = np.zeros((4, 4), dtype=object); psi[:] = Fr(0); psi[a, b] = Fr(1)
                r = ch.apply_t(k, u, psi, reg=reg)
                M[a * 4 + b] = [Fr(x) for x in r.reshape(16)]
        out.append(M)
    return out

def matvec(v, M):  # row vector times matrix
    return [sum(v[i] * M[i][j] for i in range(len(v))) for j in range(len(M[0]))]

def kron(a, b): return [x * y for x in a for y in b]

def words(v, u, ops):
    r = kron(v, u)
    x1 = matvec(r, ops[0]); x2 = matvec(r, ops[1]); x3 = matvec(r, ops[2]); x31 = matvec(x3, ops[0])
    return [r, x1, x2, x3, x31]

def rref(rows):
    """exact row reduction; returns (reduced rows, pivots)"""
    A = [list(map(Fr, r)) for r in rows]; piv = []; r0 = 0
    if not A: return A, piv
    ncol = len(A[0])
    for c in range(ncol):
        p = None
        for i in range(r0, len(A)):
            if A[i][c] != 0: p = i; break
        if p is None: continue
        A[r0], A[p] = A[p], A[r0]
        pv = A[r0][c]; A[r0] = [x / pv for x in A[r0]]
        for i in range(len(A)):
            if i != r0 and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r0])]
        piv.append(c); r0 += 1
        if r0 == len(A): break
    return A, piv

def rank(rows): return len(rref(rows)[1])

def solve_coeffs(ws, target):
    """find c with sum c_i ws_i = target (free coefficients set to 0); None if no solution"""
    n = len(ws)
    # augmented system: columns = words, rhs = target  -> solve  c . ws = target  <=> ws^T c^T = target^T
    M = [[ws[i][j] for i in range(n)] + [target[j]] for j in range(len(target))]
    R, piv = rref(M)
    if n in piv: return None
    c = [Fr(0)] * n
    for row, pc in zip(R, piv): c[pc] = row[n]
    return c

def target(v, u, kap, zz):
    """T_kappa = [v(1-kap Lambda)] (x) [u adj(1-kap Lambda)]  (polynomial normalisation)"""
    f = [1 - kap * z for z in zz]
    vv = [v[j] * f[j] for j in range(4)]
    uu = []
    for j in range(4):
        p = Fr(1)
        for m in range(4):
            if m != j: p *= f[m]
        uu.append(u[j] * p)
    return kron(vv, uu)

def kap_comps(v, u):
    """the 13 kappa-components: D = sum_a v_a e_a (x) u_a e_a, B_ab = v_a e_a (x) u_b e_b (a != b)"""
    def e(j, x):
        w = [Fr(0)] * 4; w[j] = x; return w
    D = [Fr(0)] * 16
    for j in range(4):
        D = [x + y for x, y in zip(D, kron(e(j, v[j]), e(j, u[j])))]
    comps = [D]
    for a in range(4):
        for b in range(4):
            if a != b: comps.append(kron(e(a, v[a]), e(b, u[b])))
    return comps

KAPS = [Fr(3, 10), Fr(-7, 5)]   # + 1/z4 added by callers

def test(v, u, zz, kaps=None):
    """exact test: is T_kappa in the 5-word span for each kappa?  returns (rank span, [rank(span+T)], ok, coefficients at kappa[0])"""
    ops = pair_ops_z(zz)
    pw = words(v, u, ops)
    rk = rank(pw)
    if kaps is None: kaps = KAPS + [1 / zz[3]]
    rkT = [rank(pw + [target(v, u, k, zz)]) for k in kaps]
    ok = all(x == rk for x in rkT)
    y = zz[1]
    cs = solve_coeffs(pw, [t / (1 - kaps[0] * y) ** 2 for t in target(v, u, kaps[0], zz)]) if ok else None
    return rk, rkT, ok, cs

# ---------- linear algebra helpers (exact) ----------
def nullspace(rows):
    """left-null-space helper: returns basis of {x : x . rows = 0}  (x indexed by rows)"""
    n = len(rows); m = len(rows[0])
    # solve rows^T x^T = 0  -> matrix A = rows^T (m x n)
    A = [[rows[i][j] for i in range(n)] for j in range(m)]
    R, piv = rref(A)
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        x = [Fr(0)] * n; x[f] = Fr(1)
        for row, pc in zip(R, piv): x[pc] = -row[f]
        basis.append(x)
    return basis

def inv(M):
    n = len(M); A = [list(map(Fr, M[i])) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    R, piv = rref(A)
    assert piv[:n] == list(range(n))
    return [r[n:] for r in R]

def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A): return [list(r) for r in zip(*A)]
def eye(n=4): return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
def E(i, j):
    m = [[Fr(0)] * 4 for _ in range(4)]; m[i - 1][j - 1] = Fr(1); return m
def madd(*Ms):
    return [[sum(M[i][j] for M in Ms) for j in range(len(Ms[0][0]))] for i in range(len(Ms[0]))]
def msc(c, M): return [[c * x for x in r] for r in M]

def residual_ls(pw, T, mask=None):
    """relative least-squares distance |T - c.pw|/|T| (exact rational arithmetic, float at the end); mask = list of 16 indices kept"""
    idx = list(range(len(T))) if mask is None else mask
    P_ = [[w[i] for i in idx] for w in pw]; Tm = [T[i] for i in idx]
    TT = sum(x * x for x in Tm)
    if TT == 0: return 0.0
    # independent subset of rows (pivots of transpose)
    R, piv = rref([list(c) for c in zip(*P_)])   # rows of transpose = coordinates; pivots index words
    P2 = [P_[i] for i in piv]
    G = [[sum(a * b for a, b in zip(r1, r2)) for r2 in P2] for r1 in P2]
    b = [sum(a * t for a, t in zip(r, Tm)) for r in P2]
    Gi = inv(G)
    proj = sum(b[i] * Gi[i][j] * b[j] for i in range(len(b)) for j in range(len(b)))
    res = (TT - proj) / TT
    return float(res) ** 0.5 if res > 0 else 0.0
