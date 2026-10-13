# EXACT verification of the dictionary and of the determinant overlap formula on all 16 states:
# 12 product states over Q, 4 mixed states in the number field Q[lam]/m(lam) (m = irreducible quartic)
import sympy as sp, pickle, random
U, lam = sp.symbols('U lam')
hN = sp.Rational(1, 2); th1, th2 = sp.Rational(2, 7), sp.Rational(-3, 5); th = th1 - th2
zN = [sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(4, 7)]
tops = pickle.load(open('step2.pkl', 'rb'))['tops']
cc = pickle.load(open('step1c.pkl', 'rb'))['coef']
thS, hS = sp.symbols('theta h'); vS = sp.symbols('v1:5'); uS = sp.symbols('u1:5')
V4 = lambda x: (x - th1) * (x - th2 - hN)


class K:
    """element of Q[lam]/m (or Q if m is None)"""
    m = None

    def __init__(self, p):
        P = sp.Poly(p, lam, domain='QQ')
        self.p = P.rem(K.m) if K.m is not None else P

    @staticmethod
    def c(o):
        return o if isinstance(o, K) else K(sp.sympify(o))

    def __add__(s, o): return K((s.p + K.c(o).p).as_expr())
    __radd__ = __add__

    def __sub__(s, o): return K((s.p - K.c(o).p).as_expr())

    def __rsub__(s, o): return K((K.c(o).p - s.p).as_expr())

    def __mul__(s, o): return K((s.p * K.c(o).p).as_expr())
    __rmul__ = __mul__

    def __neg__(s): return K((-s.p).as_expr())

    def inv(s):
        if K.m is None:
            return K(1 / s.p.as_expr())
        return K(sp.invert(s.p.as_expr(), K.m.as_expr(), lam))

    def __truediv__(s, o): return s * K.c(o).inv()

    def __rtruediv__(s, o): return K.c(o) * s.inv()

    def __pow__(s, n):
        r = K(1)
        for _ in range(n):
            r = r * s
        return r

    def iszero(s): return s.p.is_zero


def kdet(M):
    M = [[K.c(x) for x in r] for r in M]
    n = len(M)
    d = K(1)
    for c in range(n):
        piv = next((i for i in range(c, n) if not M[i][c].iszero()), None)
        if piv is None:
            return K(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d = d * M[c][c]
        inv = M[c][c].inv()
        for i in range(c + 1, n):
            if not M[i][c].iszero():
                f = M[i][c] * inv
                M[i] = [M[i][k] - f * M[c][k] for k in range(n)]
    return d


def ksolve_null(rows, ncols):
    M = [[K.c(x) for x in r] for r in rows]
    piv_cols = []
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(M)) if not M[i][c].iszero()), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = M[r][c].inv()
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and not M[i][c].iszero():
                f = M[i][c]
                M[i] = [M[i][k] - f * M[r][k] for k in range(ncols)]
        piv_cols.append(c)
        r += 1
    free = [c for c in range(ncols) if c not in piv_cols]
    assert len(free) == 1, ('nullity', len(free))
    fcol = free[0]
    vec = [K(0)] * ncols
    vec[fcol] = K(1)
    for i, c in enumerate(piv_cols):
        vec[c] = -M[i][fcol]
    return vec


Vv = lambda z: sp.Matrix([z**3, z**2, z, 1])


def pvec(a):
    return sp.Matrix(sp.Poly(sp.prod([U - zN[m] for m in range(4) if m != a]), U).all_coeffs())


B = sp.Matrix.hstack(*[sp.kronecker_product(pvec(i), Vv(zN[i])) for i in range(4)])
Ru = ((B.T * B).inv() * B.T * (tops[1] * B)).applyfunc(sp.expand)
C = Ru.applyfunc(lambda e: e.subs(U, 0)) - Ru[0, 0].subs(U, 0) * sp.eye(4)
mpoly = sp.Poly(C.charpoly(lam).as_expr(), lam, domain='QQ')
print('zero-weight block: char poly irreducible:', mpoly.is_irreducible, '; monic form:', sp.Poly(mpoly.as_expr() / mpoly.LC(), lam).as_expr())
adj = (C - lam * sp.eye(4)).adjugate()
xvec = adj[:, 0]
states = []
for a in range(4):
    for b in range(4):
        if a != b:
            states.append(('P%d%d' % (a + 1, b + 1), sp.kronecker_product(pvec(a), Vv(zN[b])), False))
states.append(('M(lam), all 4 mixed states at once', B * xvec, True))
random.seed(2)
results = []
for name, psi, alg in states:
    K.m = mpoly if alg else None
    psiK = [K(sp.expand(x)) for x in psi]
    n0 = 3
    assert not psiK[n0].iszero()

    def tauK(a, x):
        if a == 0:
            return K(1)
        M = tops[a].subs(U, x)
        w = K(0)
        for j in range(16):
            if M[n0, j] != 0:
                w = w + K(M[n0, j]) * psiK[j]
        return w / psiK[n0]

    Mt = tops[1].subs(U, sp.Rational(3, 11))
    lam1 = tauK(1, sp.Rational(3, 11))
    eig_ok = True
    for i in range(16):
        acc = K(0)
        for j in range(16):
            if Mt[i, j] != 0:
                acc = acc + K(Mt[i, j]) * psiK[j]
        eig_ok &= (acc - lam1 * psiK[i]).iszero()
    pts = [sp.Rational(k, 7) + sp.Rational(1, 3) for k in range(-6, 6)]
    taus_at = {x: [tauK(a, x) for a in range(5)] for x in pts}
    qs = []
    for j in range(4):
        rows = []
        for x in pts:
            row = []
            for mm in range(3):
                acc = K(0)
                for a in range(5):
                    acc = acc + taus_at[x][a] * ((-1)**a * zN[j]**(4 - a) * (x + hN - a * hN)**mm * sp.prod([V4(x - bb * hN) for bb in range(a, 4)]))
                row.append(acc)
            rows.append(row)
        qs.append(ksolve_null(rows, 3))

    def qf(j, x):
        acc = K(0)
        for mm in range(3):
            acc = acc + qs[j][mm] * x**mm
        return acc

    gA = lambda b: sp.prod([V4(th1 - l * hN) for l in range(b + 1, 4)])
    QA = lambda j, b: qf(j, th1 - b * hN) * (gA(b) * zN[j]**(-b))
    QF = lambda j, m: qf(j, th2 + m * hN) * ((V4(th2) if m == 1 else 1) * zN[j]**m)
    r = tauK(1, th2)
    ok_r = all((QF(j, 1) - r * QF(j, 0)).iszero() for j in range(4))
    ok_tau = True
    for trip in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]:
        D = [kdet([[QA(i, bb) for bb in range(4) if bb != aa] for i in trip]) for aa in range(4)]
        ok_tau &= all((D[aa] - tauK(aa, th1 - hN) * D[0]).iszero() for aa in range(4))
    ok_q1234 = kdet([[QA(i, bb) for bb in range(4)] for i in range(4)]).iszero()
    ok_vanish = all(tauK(a, th2).iszero() for a in (2, 3, 4)) and tauK(4, th1 - hN).iszero()
    ok_det = True
    ok_dec = True
    Q123 = kdet([[QA(i, bb) for bb in (1, 2, 3)] for i in range(3)])
    R = r / (hN * th)
    Rinv = R.inv()
    Qt = lambda j, bb: qf(j, th1 - bb * hN) * zN[j]**(-bb)
    Qrows = [[Qt(i, bb) for bb in range(4)] for i in range(3)]
    Qt123 = kdet([row[1:] for row in Qrows])
    cas = kdet([[Rinv**bb for bb in range(4)]] + Qrows)
    for t in range(3):
        v = [sp.Rational(random.randint(-9, 9), random.randint(1, 7)) for _ in range(4)]
        u = [sp.Rational(random.randint(-9, 9), random.randint(1, 7)) for _ in range(4)]
        sub = {thS: th, hS: hN}
        sub.update(dict(zip(vS, v)))
        sub.update(dict(zip(uS, u)))
        cnum = {k: e.subs(sub) for k, e in cc.items()}
        lhs = K(0)
        for i in range(4):
            for j in range(4):
                lhs = lhs + psiK[4 * i + j] * (v[i] * u[j])
        lhs = lhs / psiK[n0]
        w = []
        for bb in range(4):
            acc = K(0)
            for k in range(4):
                acc = acc + r**k * cnum[(bb, k)]
            w.append(acc * (-1)**bb)
        rhs = kdet([w] + [[QA(i, bb) for bb in range(4)] for i in range(3)]) / Q123
        ok_det &= (lhs - rhs).iszero()
        uR = K(0)
        for j in range(4):
            uR = uR + R**(3 - j) * u[j]
        Tt = K(0)
        for i in range(4):
            for j in range(4):
                if i < j:
                    Tt = Tt + R**(3 + i - j) * (v[i] * u[j])
        nup = []
        for bb in range(4):
            acc = K((bb * hN - th) * v[bb])
            for i in range(0, bb):
                acc = acc - Rinv**(bb - i) * (hN * v[i])
            nup.append(acc)
        rhs2 = (uR * kdet([nup] + Qrows) + Tt * cas * hN) / (Qt123 * (hN - th))
        ok_dec &= (lhs - rhs2).iszero()
    degs = [max([mm for mm in range(3) if not qs[j][mm].iszero()]) for j in range(4)]
    res = dict(state=name, eigen=eig_ok, deg_q=degs, r_all_j=ok_r, tauA_all_triples=ok_tau, Q1234_th1_zero=ok_q1234,
               t_vanish=ok_vanish, det_formula=ok_det, poly_gauge_decomp=ok_dec, R_casoratian_zero=cas.iszero())
    results.append(res)
    print(res, flush=True)
pickle.dump(results, open('exact_all.pkl', 'wb'))
