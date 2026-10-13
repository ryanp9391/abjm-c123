"""symbolic structure of the A-A matrix m_{a1 a3}(r2,r4) = sum c_{a1 k2 a3 k4} r2^k2 r4^k4 of the W_a components (needs cache/cn.pkl from run2.py)"""
import pickle, itertools, os, sympy as sp
from fractions import Fraction as Fr
WORDS = list(itertools.product(range(4), repeat=4))
cn = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'cn.pkl'), 'rb'))
r2, r4 = sp.symbols('r2 r4')
def mmat(c):
    c4 = {w: c[i] for i, w in enumerate(WORDS)}
    return sp.Matrix(4, 4, lambda a1, a3: sum(sp.Rational(c4[(a1, k2, a3, k4)].numerator, c4[(a1, k2, a3, k4)].denominator) * r2 ** k2 * r4 ** k4 for k2 in range(4) for k4 in range(4)))
for n in range(5):
    m = mmat(cn[n])
    print('=== component', n, ' rank over Q(r2,r4):', m.rank())
    # degrees
    degs = [[sp.Poly(m[i, j], r2, r4).degree_list() for j in range(4)] for i in range(4)]
    print(' degrees (r2,r4) of entries', degs)
    if n != 2:
        mins = {}
        for rows in itertools.combinations(range(4), 2):
            for cols in itertools.combinations(range(4), 2):
                d = sp.factor(m.extract(list(rows), list(cols)).det())
                mins[(rows, cols)] = d
        nz = {k: v for k, v in mins.items() if v != 0}
        print(' nonzero 2x2 minors:', len(nz), 'of 36')
        for k, v in list(nz.items())[:6]:
            print('   ', k, v)
