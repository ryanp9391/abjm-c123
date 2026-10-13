import pickle, sympy as sp
from fractions import Fraction as Fr
C, lam = pickle.load(open('phys_coeffs.pkl', 'rb'))
r2, r4 = sp.symbols('r2 r4')
# pieces: g1 = (b)-coefficient, g2 = (b^2)-coefficient in B(b) = <0|rho12(1+bE14).
# phys_coeffs has 'g^1 (a)' (A: e4->vz4 term), 'g^1 (b)', 'g^2' with gamma conventions; recombine: gamma = -b*lam1
# pair1 = (e1 - gamma vz4) x (e4 + (gamma/lam1) l1) ; terms: g^1(a) coefficient -1, g^1(b) coefficient +1, g^2 coefficient -1 (stored as +)
def poly(cdict, a1):
    return sp.expand(sum(sp.Rational(v.numerator, v.denominator) * r2**x[1] * r4**x[3] for x, v in cdict.items() if x[0] == a1 and v is not None))
for name in ['g^1 (a)', 'g^1 (b)', 'g^2']:
    print('=====', name, ' failures:', [x for x, v in C[name].items() if v is None])
    for a1 in range(4):
        P = poly(C[name], a1)
        print('  a1=%d:' % a1, sp.factor(P))
