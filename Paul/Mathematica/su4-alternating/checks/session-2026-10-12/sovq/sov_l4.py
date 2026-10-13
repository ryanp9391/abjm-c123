"""
sov_l4.py -- su(4) alternating chain A F A F (4 sites = 2 pairs), companion twist, SoV <-> Q dictionary.

CONVENTIONS (project abjm-c123, su4-alternating thread)
  * site space C^4, basis e1..e4; 256-dim chain space, flat index = 64*i1 + 16*i2 + 4*i3 + i4 (site 1 most significant,
    = Mathematica KroneckerProduct order). Bras are ROW vectors, an operator M acts on a bra as psi.M, on a ket as M.psi.
  * Lax L_s(u)[a,b] = (u-theta_s) delta_ab + h rho_s(E_ab),  rho_A(E_ab) = e_ab,  rho_F(E_ab) = -e_ba.
  * Monodromy T = L_1 L_2 L_3 L_4 (aux product, site 1 leftmost), twist on the right T^G = T G.
  * t_a(u) = sum_{|I|=a} qmin_{I,I}[T^G](u)  (rows antisymmetrised, column k at u-(k-1)h); t_4 = qdet T * det G.
    Implemented through the one-pass fused forms of dual-lax-t3.wb (checked here against the quantum minors):
      t2 = prod_s g2_s * tr[L2_1..L2_4 wedge2(G)],  g2_A = u-th-h, g2_F = u-th;  L2 = (u - th~) + h D(e_ab) rho(E_ab),
           th~_A = th, th~_F = th + h;  (wedge2 G)_{JI} = det G[J,I]
      t3 = prod_s g3_s * tr[Lb_1..Lb_4 det(G) G^{-T}],  g3_A = (u-th-h)(u-th-2h), g3_F = (u-th)(u-th-h);
           Lb_A[k,l] = (u-th+h) d_kl - h e_lk,  Lb_F[k,l] = (u-th-3h) d_kl + h e_kl.
  * Companion twist Gc(z) = {{c1,-c2,c3,-c4},{1,0,0,0},{0,1,0,0},{0,0,1,0}}, c_k = e_k(z); diagonal Lam = diag(z).
    V(z): columns (z_j^3, z_j^2, z_j, 1)^T, Gc V = V Lam.  t^{Gc}(u) R = R t^{Lam}(u), R = V^{-T} (A) (x) V (F) (x) ...
    so eigenKETS map as |Psi_c> = R |Psi_d>, bras as <x_d| = <x_c| R.
  * SoV vacuum <0| = e1 (x) e4 (x) e1 (x) e4 (companion frame).
  * SoV words <a1,k2,a3,k4| = <0| t_{a1}(th1-h) t1(th2)^k2 t_{a3}(th3-h) t1(th4)^k4.
Arithmetic: generic over Fraction (exact), float/complex (numpy) and mpmath mpf/mpc; operators are built as
dicts {(row,col): value} and converted to dense numpy arrays (object dtype for exact/mp).

SECTIONS  1-4 transfer matrices, words, rank (exact)          5 joint eigenstates (diag frame, weight blocks, 50 digits)
          6 twisted Q's (Baxter, Casoratians)                  7 SoV <-> Q dictionary (+ RV product form, N_x)
          8 ket SoV basis, Gram matrix, measure (exact)        8b dual wave functions, completeness       main()

RESULTS (verified, standard point; see sov-Q-dictionary-L4.wb for the write-up)
  * Q_i(u) = z_i^{u/h} Phi(u) q_i(u), Phi(u) = V4(u) Phi(u-h), sum_a (-1)^a tau_a(u) Q(u - a h) = 0;  deg q_I = 2 + sum_{i in I} mu_i.
  * <a1,k2,a3,k4|Psi> = <0|Psi> prod_F [Q_j(thF)/Q_j(thF-h)]^k prod_A D_a(thA)/D_0(thA),  any j, any triple I,
    D_a = det_{i in I, j in {1..4} minus {a+1}} Q_i(thA - j h);  N_x = 1 in this gauge (polynomial gauge: Nx_poly).
  * ket basis on |0bar> = e4 e1 e4 e1 (A: t1(thA)^m, F: t_n(thF + n h)): mirror dual wave functions; Gram matrix
    charge graded (S_x + S_y = 12), triangular, 1456 nonzeros -> no diagonal measure; exact measure M = G^{-1}.
"""
from fractions import Fraction as Fr
import itertools, math
import numpy as np

N = 4
TYPES = ('A', 'F', 'A', 'F')
NS = 4
DIM = N ** NS

# ---------------- standard point ----------------
H0 = Fr(1, 2)
TH0 = (Fr(2, 7), Fr(-3, 5), Fr(8, 11), Fr(-11, 8))
Z0 = (Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7))

SUBS2 = list(itertools.combinations(range(N), 2))
SUBS3 = list(itertools.combinations(range(N), 3))


def esym(z, k):
    if k == 0:
        return 1
    return sum((math.prod(c) for c in itertools.combinations(z, k)), 0 * z[0])


def Gc(z):
    c = [esym(z, k) for k in range(1, 5)]
    one, zero = 1 + 0 * z[0], 0 * z[0]
    return [[c[0], -c[1], c[2], -c[3]], [one, zero, zero, zero], [zero, one, zero, zero], [zero, zero, one, zero]]


def Gdiag(z):
    zero = 0 * z[0]
    return [[z[i] if i == j else zero for j in range(N)] for i in range(N)]


def vandermonde(z):
    return [[z[j] ** (3 - i) for j in range(N)] for i in range(N)]


# ---------- small exact/generic matrix helpers (lists of lists) ----------
def mmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum((A[i][k] * B[k][j] for k in range(m)), 0 * A[0][0]) for j in range(p)] for i in range(n)]


def mT(A):
    return [list(r) for r in zip(*A)]


def det(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    s = 0 * A[0][0]
    for j in range(n):
        if A[0][j] != 0:
            minor = [r[:j] + r[j + 1:] for r in A[1:]]
            s = s + (-1) ** j * A[0][j] * det(minor)
    return s


def minv(A):
    """inverse by Gauss-Jordan (exact for Fractions)."""
    n = len(A)
    M = [list(A[i]) + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


def wedge2(G):
    return [[det([[G[J[0]][I[0]], G[J[0]][I[1]]], [G[J[1]][I[0]], G[J[1]][I[1]]]]) for I in SUBS2] for J in SUBS2]


def cof_twist(G):
    """det G * G^{-T}"""
    d = det(G)
    Gi = minv(G)
    return [[d * Gi[j][i] for j in range(N)] for i in range(N)]


# ---------- site operators: dict {(r,c): coeff} on C^4 ----------
def rho(ty, a, b):
    return {(a, b): 1} if ty == 'A' else {(b, a): -1}


def _add(d, k, v):
    d[k] = d.get(k, 0) + v


def _wedge_coeff(k, m, I):
    if (k, m) == I:
        return 1
    if (m, k) == I:
        return -1
    return 0


# derivation action of e_ab on wedge^2:  DW[(a,b,I,J)] = D(e_ab)_{IJ}
DW = {}
for _a in range(N):
    for _b in range(N):
        _X = [[1 if (r == _a and c == _b) else 0 for c in range(N)] for r in range(N)]
        for _iI, _I in enumerate(SUBS2):
            for _iJ, _J in enumerate(SUBS2):
                _v = sum(_X[k][_J[0]] * _wedge_coeff(k, _J[1], _I) + _X[k][_J[1]] * _wedge_coeff(_J[0], k, _I) for k in range(N))
                if _v != 0:
                    DW[(_a, _b, _iI, _iJ)] = _v


def lax1(ty, u, th, h, a, b):
    d = {}
    if a == b:
        for r in range(N):
            _add(d, (r, r), u - th)
    for k, v in rho(ty, a, b).items():
        _add(d, k, h * v)
    return {k: v for k, v in d.items() if v != 0}


def lax2(ty, u, th, h, iI, iJ):
    d = {}
    if iI == iJ:
        sh = h if ty == 'F' else 0 * h
        for r in range(N):
            _add(d, (r, r), u - th - sh)
    for (a, b, I, J), c in DW.items():
        if I == iI and J == iJ:
            for k, v in rho(ty, a, b).items():
                _add(d, k, h * c * v)
    return {k: v for k, v in d.items() if v != 0}


def lax3(ty, u, th, h, k, l):
    d = {}
    if k == l:
        sh = h if ty == 'A' else -3 * h
        for r in range(N):
            _add(d, (r, r), u - th + sh)
    if ty == 'A':
        _add(d, (l, k), -h)
    else:
        _add(d, (k, l), h)
    return {kk: v for kk, v in d.items() if v != 0}


def g2(u, th, h, types=TYPES):
    return math.prod([(u - t - h) if ty == 'A' else (u - t) for ty, t in zip(types, th)])


def g3(u, th, h, types=TYPES):
    return math.prod([(u - t - h) * (u - t - 2 * h) if ty == 'A' else (u - t) * (u - t - h) for ty, t in zip(types, th)])


def transfer_dict(lax, d, twist, u, th, h, types=TYPES):
    """sum_{i,j} T_ij twist[j][i], T = L_1...L_n (aux dim d). Returns dict {(row,col): val} on the chain."""
    n = len(types)
    sitemats = [[[lax(types[s], u, th[s], h, a, b) for b in range(d)] for a in range(d)] for s in range(n)]
    P = [[dict(sitemats[0][i][k]) for k in range(d)] for i in range(d)]
    for s in range(1, n):
        Q = [[None] * d for _ in range(d)]
        for i in range(d):
            for k in range(d):
                acc = {}
                for m in range(d):
                    A = P[i][m]
                    B = sitemats[s][m][k]
                    if not A or not B:
                        continue
                    for (r1, c1), v1 in A.items():
                        for (r2, c2), v2 in B.items():
                            key = (r1 * N + r2, c1 * N + c2)
                            acc[key] = acc.get(key, 0) + v1 * v2
                Q[i][k] = acc
        P = Q
    out = {}
    for i in range(d):
        for j in range(d):
            w = twist[j][i]
            if w == 0:
                continue
            for key, v in P[i][j].items():
                out[key] = out.get(key, 0) + v * w
    return {k: v for k, v in out.items() if v != 0}


def monodromy_dicts(u, th, h, types=TYPES):
    """T_ij(u) (untwisted) as dicts, for the quantum-minor checks."""
    n = len(types)
    d = N
    sitemats = [[[lax1(types[s], u, th[s], h, a, b) for b in range(d)] for a in range(d)] for s in range(n)]
    P = [[dict(sitemats[0][i][k]) for k in range(d)] for i in range(d)]
    for s in range(1, n):
        Q = [[None] * d for _ in range(d)]
        for i in range(d):
            for k in range(d):
                acc = {}
                for m in range(d):
                    for (r1, c1), v1 in P[i][m].items():
                        for (r2, c2), v2 in sitemats[s][m][k].items():
                            key = (r1 * N + r2, c1 * N + c2)
                            acc[key] = acc.get(key, 0) + v1 * v2
                Q[i][k] = acc
        P = Q
    return P


def to_dense(dct, dtype=object, zero=0):
    if dtype is object:
        M = np.empty((DIM, DIM), dtype=object)
        M[...] = zero
    else:
        M = np.zeros((DIM, DIM), dtype=dtype)
    for (r, c), v in dct.items():
        M[r, c] = v
    return M


def identity(dtype=object, one=Fr(1), zero=Fr(0)):
    if dtype is object:
        M = np.empty((DIM, DIM), dtype=object)
        M[...] = zero
        for i in range(DIM):
            M[i, i] = one
        return M
    return np.identity(DIM, dtype=dtype)


class Chain:
    """Transfer matrices of the 4-site A F A F chain for given (h, theta, twist G).
    exact=True: Fraction entries (object arrays); exact=False: complex128 numpy."""

    def __init__(self, h=H0, th=TH0, G=None, z=Z0, types=TYPES, exact=True):
        self.h, self.th, self.types = h, tuple(th), types
        self.G = G if G is not None else Gc(z)
        self.G2 = wedge2(self.G)
        self.G3 = cof_twist(self.G)
        self.detG = det(self.G)
        self.exact = exact
        self.dtype = object if exact else complex
        self._cache = {}

    def t(self, a, u, reduced=False):
        key = (a, u, reduced)
        if key in self._cache:
            return self._cache[key]
        h, th = self.h, self.th
        if a == 0:
            M = identity(self.dtype)
        elif a == 1:
            M = to_dense(transfer_dict(lax1, 4, self.G, u, th, h, self.types), self.dtype)
        elif a == 2:
            M = to_dense(transfer_dict(lax2, 6, self.G2, u, th, h, self.types), self.dtype)
            if not reduced:
                M = M * g2(u, th, h, self.types)
        elif a == 3:
            M = to_dense(transfer_dict(lax3, 4, self.G3, u, th, h, self.types), self.dtype)
            if not reduced:
                M = M * g3(u, th, h, self.types)
        else:
            raise ValueError(a)
        self._cache[key] = M
        return M


def product_vec(vs):
    """tensor product of 4 site vectors -> flat object array"""
    out = [1]
    for v in vs:
        out = [a * b for a in out for b in v]
    return np.array(out, dtype=object)


def unit(i, zero=Fr(0), one=Fr(1)):
    return [one if k == i else zero for k in range(N)]


VAC_SITES = (0, 3, 0, 3)  # A: e1, F: e4 (0-based)


def sov_vacuum():
    return product_vec([unit(i) for i in VAC_SITES])


WORDS = list(itertools.product(range(4), repeat=4))  # (a1,k2,a3,k4)


def sov_words(ch, vac=None, pts=None):
    """rows <a1,k2,a3,k4| = <0| t_{a1}(p1) t1(p2)^k2 t_{a3}(p3) t1(p4)^k4, default p = (th1-h, th2, th3-h, th4)."""
    h, th = ch.h, ch.th
    if pts is None:
        pts = (th[0] - h, th[1], th[2] - h, th[3])
    vac = sov_vacuum() if vac is None else vac
    A1 = [vac.dot(ch.t(a, pts[0])) if a else vac for a in range(4)]
    T2 = ch.t(1, pts[1])
    T3 = [ch.t(a, pts[2]) for a in range(4)]
    T4 = ch.t(1, pts[3])
    rows = []
    cache = {}
    for (a1, k2, a3, k4) in WORDS:
        if (a1, k2) not in cache:
            v = A1[a1]
            for _ in range(k2):
                v = v.dot(T2)
            cache[(a1, k2)] = v
        v = cache[(a1, k2)]
        v = v.dot(T3[a3]) if a3 else v
        for _ in range(k4):
            v = v.dot(T4)
        rows.append(v)
    return np.array(rows)


def rank_mod_p(M, p=2147483629):
    """exact rank of a rational matrix (object array of Fractions) modulo a large prime (lower bound of the rational rank;
    equal unless p divides a minor)."""
    A = np.array([[(x.numerator % p) * pow(x.denominator % p, -1, p) % p for x in row] for row in M], dtype=object)
    A = A.astype(np.int64)
    n, m = A.shape
    r = 0
    for c in range(m):
        piv = None
        for i in range(r, n):
            if A[i, c] % p != 0:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), -1, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nz = np.nonzero(col)[0]
        for i in nz:
            A[i] = (A[i] - int(col[i]) * A[r]) % p
        r += 1
        if r == n:
            break
    return r


# =====================================================================================================================
# 5. Joint eigenstates (diagonal frame, weight blocks, mpmath) and eigenvalue polynomials
# =====================================================================================================================
import mpmath as mp
import pickle, os, time

mp.mp.dps = 50


def tomp(x):
    """Fraction/int/float -> mpf"""
    if isinstance(x, Fr):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)


def weight(idx):
    """Cartan weight of a basis state (i1,j2,i3,j4): A site +e_i, F site -e_j (vacuum e1 e4 e1 e4 -> (2,0,0,-2))"""
    w = [0] * 4
    for s, i in enumerate(idx):
        w[i] += 1 if TYPES[s] == 'A' else -1
    return tuple(w)


def weight_blocks():
    blocks = {}
    for n, idx in enumerate(itertools.product(range(4), repeat=4)):
        blocks.setdefault(weight(idx), []).append(n)
    return blocks


UPTS = [Fr(k, 3) + Fr(1, 17) for k in range(-2, 3)]          # interpolation points for the reduced tau-bar_a (deg 4)
COMB = ((1, Fr(3, 11), 1), (3, Fr(-5, 7), Fr(7, 10)), (2, Fr(2, 9), Fr(-3, 13)))   # generic combination to diagonalise


def eigenstates(dps=50):
    """All 256 joint eigenstates of t1, t2, t3 with diagonal twist Lam = diag(z) (equivalently of t^{Gc}, |Psi_c> = R|Psi_d>).
    Returns list of dicts: weight, idx (block basis indices), vec (mp, block components), taubar {a: coeffs low->high},
    res (max eigen-residual)."""
    mp.mp.dps = dps
    chd = Chain(G=Gdiag(Z0))
    mats = {(a, u): chd.t(a, u, reduced=True) for a in (1, 2, 3) for u in UPTS}
    Mc = COMB[0][2] * chd.t(1, COMB[0][1])
    for a, u, c in COMB[1:]:
        Mc = Mc + c * chd.t(a, u, reduced=True)
    blocks = weight_blocks()
    Vm = mp.matrix([[tomp(u) ** j for j in range(5)] for u in UPTS])
    out = []
    for w, idx in sorted(blocks.items()):
        B = mp.matrix([[tomp(Mc[r, c]) for c in idx] for r in idx])
        E, ER = mp.eig(B)
        bm = {k: mp.matrix([[tomp(M[r, c]) for c in idx] for r in idx]) for k, M in mats.items()}
        for k in range(len(idx)):
            v = ER[:, k]
            p = max(range(len(idx)), key=lambda r: abs(v[r]))
            v = v / v[p]
            tb, res = {}, mp.mpf(0)
            for a in (1, 2, 3):
                vals = []
                for u in UPTS:
                    tv = bm[(a, u)] * v
                    lam = tv[p] / v[p]
                    res = max(res, mp.norm(tv - lam * v))
                    vals.append(lam)
                c = mp.lu_solve(Vm, mp.matrix(vals))
                tb[a] = [c[j] for j in range(5)]
            out.append({'weight': w, 'idx': idx, 'vec': [v[r] for r in range(len(idx))], 'taubar': tb, 'res': res,
                        'combEig': E[k]})
    return out


# vacuum data (generic theta per site), mp
_h, _th, _z = tomp(H0), [tomp(t) for t in TH0], [tomp(t) for t in Z0]


def Vvac(i, u, h=None, th=None):
    """eigenvalue of T_ii(u) (i=1..4) on the ket vacuum e1 e4 e1 e4 (diag frame):
    V1 = prod_A (u-thA+h) prod_F (u-thF), V2 = V3 = prod_A (u-thA) prod_F (u-thF), V4 = prod_A (u-thA) prod_F (u-thF-h)"""
    h = _h if h is None else h
    th = _th if th is None else th
    A = [u - th[s] for s in range(4) if TYPES[s] == 'A']
    F = [u - th[s] for s in range(4) if TYPES[s] == 'F']
    if i == 1:
        return mp.fprod([x + h for x in A]) * mp.fprod(F)
    if i in (2, 3):
        return mp.fprod(A) * mp.fprod(F)
    return mp.fprod(A) * mp.fprod([x - h for x in F])


def g2mp(u):
    return mp.fprod([(u - t - _h) if ty == 'A' else (u - t) for ty, t in zip(TYPES, _th)])


def g3mp(u):
    return mp.fprod([(u - t - _h) * (u - t - 2 * _h) if ty == 'A' else (u - t) * (u - t - _h) for ty, t in zip(TYPES, _th)])


def qdet_site(ty, u, t):
    """one-site quantum determinant: A (u-t+h)(u-t-h)(u-t-2h)(u-t-3h), F (u-t)(u-t-h)(u-t-2h)(u-t-4h)"""
    h = _h
    if ty == 'A':
        return (u - t + h) * (u - t - h) * (u - t - 2 * h) * (u - t - 3 * h)
    return (u - t) * (u - t - h) * (u - t - 2 * h) * (u - t - 4 * h)


def peval(c, u):
    r = mp.mpf(0)
    for x in reversed(c):
        r = r * u + x
    return r


class Eig:
    """one joint eigenstate: full eigenvalues tau_a(u) (a = 0..4), tau_4 = det(G) qdet"""

    def __init__(self, d):
        self.d, self.weight, self.tb = d, d['weight'], d['taubar']

    def tau(self, a, u):
        if a == 0:
            return mp.mpf(1)
        if a == 1:
            return peval(self.tb[1], u)
        if a == 2:
            return g2mp(u) * peval(self.tb[2], u)
        if a == 3:
            return g3mp(u) * peval(self.tb[3], u)
        return mp.fprod(_z) * mp.fprod([qdet_site(ty, u, t) for ty, t in zip(TYPES, _th)])


# =====================================================================================================================
# 6. Twisted Q-functions
#    Baxter: sum_{a=0}^4 (-1)^a tau_a(u) Q(u - a h) = 0  (this module's argument convention; RV/coordinator: Q^RV(u) = Q(u-h)).
#    On the vacuum the Baxter operator is (1 - z1V1 D)(1 - z2V2 D)(1 - z3V3 D)(1 - z4V4 D), D f(u) = f(u-h).
#    Gauge: Q_i(u) = z_i^{u/h} Phi(u) q_i(u), Phi(u) = V4(u) Phi(u-h) (RV's Gamma[nu_1] gauge), q_i monic polynomial:
#        sum_a (-1)^a z_i^{-a} tau_a(u) q_i(u - a h) prod_{b=a}^{3} V4(u - b h) = 0.
#    deg q_i = 2 + mu_i (mu = weight).  Casoratians Q_I(u) = det_{i in I, s<|I|} Q_i(u - s h); polynomial parts
#    P_I(u) = det[z_i^{|I|-1-s} q_i(u - s h)] = F_|I|(u) q_I(u),  F_2 = prod_F (u-thF), F_3 = prod_F (u-thF)(u-thF-h),
#    F_4 = prod_F (u-thF)(u-thF-h)(u-thF-2h) prod_A (u-thA+h);  q_1234 = const (state independent).
# =====================================================================================================================
def _padd(a, b):
    n = max(len(a), len(b))
    return [(a[k] if k < len(a) else 0) + (b[k] if k < len(b) else 0) for k in range(n)]


def _pmul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def _pshift(a, c):
    out = [mp.mpf(0)] * len(a)
    for k, ak in enumerate(a):
        for j in range(k + 1):
            out[j] += ak * mp.binomial(k, j) * c ** (k - j)
    return out


def _ptrim(a, tol=mp.mpf(10) ** -30):
    s = max(abs(x) for x in a)
    b = list(a)
    while len(b) > 1 and abs(b[-1]) <= tol * s:
        b.pop()
    return b


def _pdet(M):
    n = len(M)
    tot = [mp.mpf(0)]
    for p in itertools.permutations(range(n)):
        s = 1
        for i in range(n):
            for j in range(i + 1, n):
                if p[i] > p[j]:
                    s = -s
        term = [mp.mpf(s)]
        for i in range(n):
            term = _pmul(term, M[i][p[i]])
        tot = _padd(tot, term)
    return tot


def _pdivmod(a, b):
    a = list(a)
    q = [mp.mpf(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b):
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        q[k] = c
        for j in range(len(b)):
            a[k + j] -= c * b[j]
        a.pop()
    return q, a


def _from_roots(rs):
    p = [mp.mpf(1)]
    for r in rs:
        p = _pmul(p, [-r, mp.mpf(1)])
    return p


_SAMPLE = [mp.mpf(k) / 7 + mp.mpf(1) / 29 for k in range(-20, 21)]


def solve_q(st, i, Dmax=10, tol=mp.mpf(10) ** -35):
    """minimal-degree monic polynomial q_i (coefficients low->high) of the twisted Baxter equation for twist z_i.
    Rows are scaled by the size of the individual Baxter terms (a small singular value = genuine cancellation)."""
    zi = _z[i - 1]
    for D in range(Dmax + 1):
        rows = []
        for u in _SAMPLE[:22 + D]:
            row, scale = [mp.mpf(0)] * (D + 1), mp.mpf(0)
            for a in range(5):
                pre = (-1) ** a * zi ** (-a) * st.tau(a, u) * mp.fprod([Vvac(4, u - b * _h) for b in range(a, 4)])
                x = u - a * _h
                for k in range(D + 1):
                    row[k] += pre * x ** k
                scale = max(scale, abs(pre) * max(abs(x) ** k for k in range(D + 1)))
            rows.append([r / scale for r in row])
        U, S, Vt = mp.svd_c(mp.matrix(rows))
        if S[S.rows - 1] < tol:
            v = [mp.conj(Vt[Vt.rows - 1, c]) for c in range(Vt.cols)]
            return [x / v[-1] for x in v], S[S.rows - 1]
    raise RuntimeError('no polynomial solution')


THA = [_th[s] for s in range(4) if TYPES[s] == 'A']
THF = [_th[s] for s in range(4) if TYPES[s] == 'F']
FVAC = {2: _from_roots(THF), 3: _from_roots(THF + [t + _h for t in THF]),
        4: _from_roots(THF + [t + _h for t in THF] + [t + 2 * _h for t in THF] + [t - _h for t in THA])}
SUBSETS = [I for k in (1, 2, 3, 4) for I in itertools.combinations((1, 2, 3, 4), k)]


def q_functions(st):
    """all 15 q_I (monic singles; Casoratian polynomial parts divided by the vacuum factors F_|I|).
    Returns (dict I -> coeffs, max relative remainder of the divisions, max Baxter singular value)"""
    Q, rem, sv = {}, mp.mpf(0), mp.mpf(0)
    for i in (1, 2, 3, 4):
        Q[(i,)], s = solve_q(st, i)
        sv = max(sv, s)
    for I in SUBSETS:
        if len(I) == 1:
            continue
        k = len(I)
        M = [[[_z[m - 1] ** (k - 1 - s) * c for c in _pshift(Q[(m,)], -s * _h)] for m in I] for s in range(k)]
        P = _ptrim(_pdet(M))
        q, r = _pdivmod(P, FVAC[k])
        rem = max(rem, max([abs(x) for x in r] + [0]) / max(abs(x) for x in P))
        Q[I] = _ptrim(q)
    return Q, rem, sv


# =====================================================================================================================
# 7. The SoV <-> Q dictionary (verified for all 256 words x 256 states)
#    Regularised Gamma-gauge values at site s (word-independent factor z_i^{th_s/h} Phi(th_s - 5h) dropped):
#        Qt_i(th_s - j h) = z_i^{-j} prod_{b=j}^{4} V4(th_s - b h) q_i(th_s - j h).
#    F site (t1 powers at thF):  tau_1(thF) = Q_j(thF)/Q_j(thF - h)  for ANY j = 1..4           (N = 1)
#    A site (t_a at thA - h):    tau_a(thA - h) = D_a/D_0,  D_a = det_{i in I; j in {1,2,3,4}\{a+1}} Q_i(thA - j h)
#                                for ANY triple I (D_0 = Q_I(thA - 2h), D_3 = Q_I(thA - h): Casoratians)  (N = 1)
#    => <a1,k2,a3,k4|Psi> = <0|Psi> * prod_F [Q(thF)/Q(thF-h)]^k * prod_A D_a/D_0   (no renormalisation needed in this gauge)
#    Polynomial gauge (qh_i(u) = z_i^{u/h} q_i(u)): <x|Psi>/<0|Psi> = N_x * prod_F [qh(thF)/qh(thF-h)]^k
#        * prod_A det[qh_i(thA - j h)]_{j in S_a} / det[qh_i(thA - j h)]_{j in {2,3,4}},
#        N_x = prod_F V4(thF)^k * prod_A prod_{b=1}^{a} V4(thA - b h)   (exact rationals, see Nx_poly).
#    RV form: with x the B-operator zeros (F: x in {thF, thF+h}^3, k of them at thF+h; A: {thA..thA-3h} minus thA - a h,
#    plus frozen ones) <Lambda|Psi>/<0|Psi> = prod_s prod_{k<=3} det_{i,j<=k} Q^RV_{sigma(i)}(x_kj) / (same at Lambda = 0),
#    Q^RV(u) = Q(u - h), any sigma.
# =====================================================================================================================
def Qt(Q, i, ths, j, J=5):
    return _z[i - 1] ** (-j) * mp.fprod([Vvac(4, ths - b * _h) for b in range(j, J)]) * peval(Q[(i,)], ths - j * _h)


def wave_site(Q, s, lab, j=1, I=(1, 2, 3)):
    """site factor of <x|Psi>/<0|Psi>: F: (Q_j(th)/Q_j(th-h))^k; A: D_a/D_0 with the triple I"""
    ths = _th[s]
    if TYPES[s] == 'F':
        return (Qt(Q, j, ths, 0) / Qt(Q, j, ths, 1)) ** lab
    D = lambda pts: mp.det(mp.matrix([[Qt(Q, i, ths, jj) for jj in pts] for i in I]))
    return D([jj for jj in (1, 2, 3, 4) if jj != lab + 1]) / D((2, 3, 4))


def wave(Q, word, j=1, I=(1, 2, 3)):
    return mp.fprod([wave_site(Q, s, word[s], j, I) for s in range(4)])


def tau_word(st, word):
    """<word|Psi>/<0|Psi> = tau_{a1}(th1-h) tau1(th2)^k2 tau_{a3}(th3-h) tau1(th4)^k4"""
    return (st.tau(word[0], _th[0] - _h) * st.tau(1, _th[1]) ** word[1] * st.tau(word[2], _th[2] - _h)
            * st.tau(1, _th[3]) ** word[3])


def _rv_points(ty, lab):
    """RV separated variables in units of j (this module's Q argument = th - j h), levels k = 1,2,3"""
    if ty == 'F':
        mu = [1] * lab + [0] * (3 - lab)
        xs = [[mu[0]], [mu[1], -1], [mu[2], -1, -2]]
    else:
        xs = [[0], [0, -1], [v for v in (0, -1, -2, -3) if v != -lab]]
    return [[1 - v for v in lev] for lev in xs]


def rv_site_table(Q, sigma=(1, 2, 3)):
    """[s][lab] -> RV2 (5.20) site factor prod_k det_{i,j<=k} Q_sigma(i)(x_kj) / (same at lab = 0), regularised gauge"""
    def site(s, lab):
        p = mp.mpf(1)
        for k, js in enumerate(_rv_points(TYPES[s], lab), start=1):
            p *= mp.det(mp.matrix([[Qt(Q, sigma[i], _th[s], jj) for jj in js] for i in range(k)]))
        return p
    tab = []
    for s in range(4):
        vals = [site(s, lab) for lab in range(4)]
        tab.append([v / vals[0] for v in vals])
    return tab


def rv_product(Q, word, sigma=(1, 2, 3), table=None):
    """RV2 (5.20) product prod_s prod_k det_{i,j<=k} Q_sigma(i)(x_kj), divided by its vacuum value"""
    tab = rv_site_table(Q, sigma) if table is None else table
    return mp.fprod([tab[s][word[s]] for s in range(4)])


def V4_exact(u, h=H0, th=TH0):
    return math.prod([(u - t) if ty == 'A' else (u - t - h) for ty, t in zip(TYPES, th)])


def Nx_poly(word, h=H0, th=TH0):
    """exact normalisation of the words in the polynomial gauge (see header of section 7)"""
    n = Fr(1)
    for s, lab in enumerate(word):
        if TYPES[s] == 'F':
            n *= V4_exact(th[s]) ** lab
        else:
            n *= math.prod([V4_exact(th[s] - b * h) for b in range(1, lab + 1)])
    return n


# =====================================================================================================================
# 8. Ket SoV basis, dual basis, Gram matrix and measure
#    Ket vacuum |0bar> = e4 (x) e1 (x) e4 (x) e1;
#    ket words |m1,n2,m3,n4> = t1(th1)^m1 t_{n2}(th2 + n2 h) t1(th3)^m3 t_{n4}(th4 + n4 h)|0bar>.
#    Dual wave functions (dual Q's w_I(u) = Z_I^{(u-th)/h} q_I(u), Z_I = prod_{i in I} z_i, any triple(s)):
#      A: tau1(thA) = NA * w_I(thA - 2h)/w_I(thA - h);
#      F: tau_n(thF + n h) = NF_n det[w_I(thF+(j-1)h)]_{j != n} / det[...]_{j=1,2,3}.
#    Gram G_xy = <x|y>: nonzero only if sum(x) + sum(y) = 12, always nonzero for y = 3 - x, triangular (acyclic),
#    1456 nonzeros; measure M = G^{-1} (same pattern): 1 = sum_{x,y} |y> M_{yx} <x|.
# =====================================================================================================================
KET_VAC = (3, 0, 3, 0)


def ket_words(ch):
    """columns |y> (object array 256x256), y in WORDS order"""
    vac = product_vec([unit(i) for i in KET_VAC])
    h, th = ch.h, ch.th
    cols = []
    for lab in WORDS:
        v = vac
        for s in (3, 2, 1, 0):
            if TYPES[s] == 'A':
                for _ in range(lab[s]):
                    v = ch.t(1, th[s]).dot(v)
            elif lab[s]:
                v = ch.t(lab[s], th[s] + lab[s] * h).dot(v)
        cols.append(v)
    return np.array(cols).T


def gram(W, Y):
    """exact G[i, j] = <word_i | ket_j>"""
    Wnz = [[(c, x) for c, x in enumerate(r) if x != 0] for r in W]
    G = np.empty((DIM, DIM), dtype=object)
    for i in range(DIM):
        for j in range(DIM):
            G[i, j] = sum((x * Y[c, j] for c, x in Wnz[i]), Fr(0))
    return G


def measure(G):
    """exact M = G^{-1} by charge blocks. Returns dict {(y, x): M_yx} (labels as WORDS tuples, y = ket label)."""
    idx = {w: i for i, w in enumerate(WORDS)}
    blocks = {}
    for x in WORDS:
        blocks.setdefault(sum(x), []).append(x)
    M = {}
    for S, xs in blocks.items():
        ys = [tuple(3 - c for c in x) for x in xs]          # kets with charge 12 - S
        B = [[G[idx[x], idx[y]] for y in ys] for x in xs]
        Bi = minv(B)                                         # rows: kets y, cols: bras x
        for i, y in enumerate(ys):
            for j, x in enumerate(xs):
                if Bi[i][j] != 0:
                    M[(y, x)] = Bi[i][j]
    return M


def is_triangular(G):
    """True if the x <-> 3-x relabelled Gram matrix is triangular w.r.t. some order (off-diagonal graph acyclic)"""
    idx = {w: i for i, w in enumerate(WORDS)}
    adj = {x: [xp for xp in WORDS if xp != x and G[idx[x], idx[tuple(3 - c for c in xp)]] != 0] for x in WORDS}
    state = {}

    def dfs(u):
        state[u] = 1
        for v in adj[u]:
            if state.get(v) == 1 or (state.get(v) is None and dfs(v)):
                return True
        state[u] = 2
        return False
    return not any(dfs(u) for u in WORDS if state.get(u) is None)


# =====================================================================================================================
# 8b. Dual wave functions and completeness on the eigenstates
# =====================================================================================================================
def dual_w(Q, I, u, ths):
    """dual (triple-index) Q value w_I(u) = Z_I^{(u-ths)/h} q_I(u)"""
    return mp.fprod([_z[i - 1] for i in I]) ** ((u - ths) / _h) * peval(Q[I], u)


TRIPLES = list(itertools.combinations((1, 2, 3, 4), 3))


def dual_constants(st, Q):
    """A sites: NA_s = tau1(th_s) w_I(th_s-h)/w_I(th_s-2h);  F sites: NF_{s,n} = tau_n(th_s+nh) det_{j=1,2,3}/det_{j!=n}
    (w at th_s + (j-1)h, three triples).  State independent (checked in main)."""
    out = {}
    for s in range(4):
        if TYPES[s] == 'A':
            out[(s, 1)] = [st.tau(1, _th[s]) * dual_w(Q, I, _th[s] - _h, _th[s]) / dual_w(Q, I, _th[s] - 2 * _h, _th[s]) for I in TRIPLES]
        else:
            mm = [[dual_w(Q, I, _th[s] + (j - 1) * _h, _th[s]) for j in range(4)] for I in TRIPLES[:3]]
            for a in (1, 2, 3):
                Dn = mp.det(mp.matrix([[mm[l][j] for j in range(4) if j != a] for l in range(3)]))
                D0 = mp.det(mp.matrix([[mm[l][j] for j in (1, 2, 3)] for l in range(3)]))
                out[(s, a)] = [st.tau(a, _th[s] + a * _h) * D0 / Dn]
    return out


def _apply_kron(mats, vec, left=False):
    v = list(vec)
    for s in range(4):
        m = mats[s]
        new = [mp.mpf(0)] * DIM
        stride = 4 ** (3 - s)
        for idx in range(DIM):
            if v[idx] == 0:
                continue
            d = (idx // stride) % 4
            base = idx - d * stride
            for k in range(4):
                c = m[d][k] if left else m[k][d]
                if c != 0:
                    new[base + k * stride] += tomp(c) * v[idx]
        v = new
    return v


def companion_vectors(eig):
    """right eigenkets |Psi_c> = R|Psi_d> and left eigenbras <PsiL_c| = <PsiL_d| R^{-1} (normalised <PsiL|Psi> = 1)"""
    V = vandermonde(Z0)
    VmT, Vi, VT = mT(minv(V)), minv(V), mT(V)
    byblock = {}
    for n, d in enumerate(eig):
        byblock.setdefault(tuple(d['idx']), []).append(n)
    right, left = [None] * len(eig), [None] * len(eig)
    for idx, ns in byblock.items():
        Rb = mp.matrix([[eig[n]['vec'][r] for n in ns] for r in range(len(idx))])
        Lb = mp.inverse(Rb)
        for k, n in enumerate(ns):
            v = [mp.mpf(0)] * DIM
            l = [mp.mpf(0)] * DIM
            for r, c in enumerate(idx):
                v[c] = Rb[r, k]
                l[c] = Lb[k, r]
            right[n] = _apply_kron([VmT, V, VmT, V], v)
            left[n] = _apply_kron([VT, Vi, VT, Vi], l, left=True)
    return right, left


# =====================================================================================================================
# main: reproduce the checks (python half of sov-Q-dictionary-L4.wb). Caches in ./cache_sov_l4/
# =====================================================================================================================
def main(cache='cache_sov_l4'):
    os.makedirs(cache, exist_ok=True)
    T0 = time.time()
    rep = {}

    def cached(name, fn):
        f = os.path.join(cache, name + '.pkl')
        if os.path.exists(f):
            return pickle.load(open(f, 'rb'))
        r = fn()
        pickle.dump(r, open(f, 'wb'))
        return r
    ch = Chain()
    W = cached('words', lambda: sov_words(ch))
    rep['rank words (mod p, exact)'] = rank_mod_p(W)
    eig = cached('eigen', eigenstates)
    sts = [Eig(d) for d in eig]
    rep['# eigenstates'] = len(sts)
    rep['max eigen residual'] = float(max(d['res'] for d in eig))
    QQ = cached('Q', lambda: [q_functions(st) for st in sts])
    rep['max Casoratian remainder'] = float(max(q[1] for q in QQ))
    rep['max Baxter singular value'] = float(max(q[2] for q in QQ))
    rep['deg q_i = 2 + mu_i for all states'] = all(len(q[0][(i,)]) - 1 == 2 + st.weight[i - 1]
                                                   for q, st in zip(QQ, sts) for i in (1, 2, 3, 4))
    q1234 = [q[0][(1, 2, 3, 4)][0] for q in QQ]
    rep['q_1234 spread'] = float(max(abs(x / q1234[0] - 1) for x in q1234))
    errF = errA = errRV = errN = 0
    for st, (Q, _, _) in zip(sts, QQ):
        for j in (1, 2, 3, 4):
            for s in (1, 3):
                errF = max(errF, abs(wave_site(Q, s, 1, j=j) / st.tau(1, _th[s]) - 1))
        for I in itertools.combinations((1, 2, 3, 4), 3):
            for s in (0, 2):
                for a in (1, 2, 3):
                    errA = max(errA, abs(wave_site(Q, s, a, I=I) / st.tau(a, _th[s] - _h) - 1))
        for sigma in [(1, 2, 3), (4, 3, 2), (2, 4, 1), (3, 1, 4)]:
            tab = rv_site_table(Q, sigma)
            for w in WORDS:
                errRV = max(errRV, abs(rv_product(Q, w, table=tab) / tau_word(st, w) - 1))
        ntab = []
        for s in range(4):
            row = []
            for lab in range(4):
                if TYPES[s] == 'F':
                    row.append((tomp(V4_exact(TH0[s])) * _z[3] * peval(Q[(4,)], _th[s]) / peval(Q[(4,)], _th[s] - _h)) ** lab)
                else:
                    qh = lambda i, jj: _z[i - 1] ** (-jj) * peval(Q[(i,)], _th[s] - jj * _h)
                    Sa = [jj for jj in (1, 2, 3, 4) if jj != lab + 1]
                    num = mp.det(mp.matrix([[qh(i, jj) for jj in Sa] for i in (1, 2, 3)]))
                    den = mp.det(mp.matrix([[qh(i, jj) for jj in (2, 3, 4)] for i in (1, 2, 3)]))
                    row.append(tomp(math.prod([V4_exact(TH0[s] - b * H0) for b in range(1, lab + 1)])) * num / den)
            ntab.append(row)
        for w in WORDS:
            val = mp.fprod([ntab[s][w[s]] for s in range(4)])
            errN = max(errN, abs(val / tau_word(st, w) - 1))
    rep['F site: tau1(thF) = Q_j(thF)/Q_j(thF-h), all j'] = float(errF)
    rep['A site: tau_a(thA-h) = D_a/D_0, all triples'] = float(errA)
    rep['RV product formula, all words, 4 sigma'] = float(errRV)
    rep['polynomial gauge with exact N_x, all words'] = float(errN)
    Y = cached('kets', lambda: ket_words(ch))
    rep['rank ket words (mod p)'] = rank_mod_p(Y.T)
    G = cached('gram', lambda: gram(W, Y))
    nz = [(WORDS[i], WORDS[j]) for i in range(DIM) for j in range(DIM) if G[i, j] != 0]
    rep['Gram nonzeros'] = len(nz)
    rep['Gram charge rule S_x + S_y = 12'] = all(sum(x) + sum(y) == 12 for x, y in nz)
    rep['Gram triangular'] = is_triangular(G)
    M = cached('measure', lambda: measure(G))
    rep['measure nonzeros'] = len(M)
    # dual wave functions: constants state independent
    consts = [dual_constants(st, Q) for st, (Q, _, _) in zip(sts, QQ)]
    spread = 0
    for key in consts[0]:
        ref = consts[0][key][0]
        for c in consts:
            spread = max(spread, max(abs(x / ref - 1) for x in c[key]))
    rep['dual wave functions: spread of N_A, N_F over states/triples'] = float(spread)
    # completeness with the exact measure on a sample of eigenstates
    right, left = cached('lr', lambda: companion_vectors(eig))
    Wm = [[(c, tomp(x)) for c, x in enumerate(r) if x != 0] for r in W]
    Ym = [[(c, tomp(Y[c, j])) for c in range(DIM) if Y[c, j] != 0] for j in range(DIM)]
    idx = {w: i for i, w in enumerate(WORDS)}
    errC = 0
    for n in range(0, 256, 17):
        ov = [mp.fsum(x * right[n][c] for c, x in Wm[i]) for i in range(DIM)]
        ovL = [mp.fsum(left[n][c] * x for c, x in Ym[j]) for j in range(DIM)]
        tot = mp.fsum(ovL[idx[y]] * tomp(v) * ov[idx[x]] for (y, x), v in M.items())
        nrm = mp.fsum(left[n][c] * right[n][c] for c in range(DIM))
        errC = max(errC, abs(tot / nrm - 1))
    rep['completeness 1 = sum |y> M_yx <x| on 16 eigenstates'] = float(errC)
    rep['time (s)'] = time.time() - T0
    for k, v in rep.items():
        print('%-55s %s' % (k, v))
    return rep


if __name__ == '__main__':
    main()
