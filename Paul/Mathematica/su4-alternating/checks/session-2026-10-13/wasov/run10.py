import pickle, itertools, os, functools, sympy as sp
gcd_list = lambda L: functools.reduce(lambda a, b: sp.gcd(a, b), L)
WORDS = list(itertools.product(range(4), repeat=4))
cn = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'cn.pkl'), 'rb'))
r2, r4 = sp.symbols('r2 r4')
def mmat(c):
    c4 = {w: c[i] for i, w in enumerate(WORDS)}
    return sp.Matrix(4, 4, lambda a1, a3: sum(sp.Rational(c4[(a1, k2, a3, k4)].numerator, c4[(a1, k2, a3, k4)].denominator) * r2 ** k2 * r4 ** k4 for k2 in range(4) for k4 in range(4)))
def bideg(p): 
    P = sp.Poly(p, r2, r4); return P.degree(r2), P.degree(r4)
for n in (0, 1, 3, 4):
    m = mmat(cn[n])
    g1 = gcd_list([sp.Poly(m[i, j], r2, r4) for i in range(4) for j in range(4)])
    mins = [m.extract(list(R), list(C)).det() for R in itertools.combinations(range(4), 2) for C in itertools.combinations(range(4), 2)]
    g2 = gcd_list([sp.Poly(sp.expand(x), r2, r4) for x in mins])
    # row-Plucker gcd (common factor over columns for fixed rows) and column-Plucker gcd
    fR = {}
    for R in itertools.combinations(range(4), 2):
        fR[R] = gcd_list([sp.Poly(sp.expand(m.extract(list(R), list(C)).det()), r2, r4) for C in itertools.combinations(range(4), 2)])
    gC = {}
    for C in itertools.combinations(range(4), 2):
        gC[C] = gcd_list([sp.Poly(sp.expand(m.extract(list(R), list(C)).det()), r2, r4) for R in itertools.combinations(range(4), 2)])
    print('component', n, ': gcd of entries bidegree', bideg(g1.as_expr()) if g1.degree() > 0 else 'const',
          '| gcd of 2x2 minors bidegree', bideg(g2.as_expr()),
          '| row-Pluecker gcds bidegrees', [bideg(v.as_expr()) if v.degree() > 0 else None for v in fR.values()],
          '| col-Pluecker gcds', [bideg(v.as_expr()) if v.degree() > 0 else None for v in gC.values()], flush=True)
    for R, v in list(fR.items())[:2]:
        print('   fR', R, sp.factor_list(v.as_expr())[1][:3] if v.degree() > 0 else 1)
