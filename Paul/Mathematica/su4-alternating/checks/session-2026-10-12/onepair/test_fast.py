# sanity: fast code reproduces the 2-site SoV matrix W (mod p) from the symbolic one
import sympy as sp, pickle, numpy as np
from fractions import Fraction
from fastchain import *
th1, th2, h = sp.symbols('theta1 theta2 h')
d = pickle.load(open('step1.pkl', 'rb')); W = d['W']
hN = Fraction(1, 2); thN = [Fraction(2, 7), Fraction(-3, 5)]
zN = [Fraction(3, 2), Fraction(5, 4), Fraction(2, 3), Fraction(4, 7)]
from itertools import combinations
cvec = [sum(np.prod(s) for s in combinations(zN, k)) for k in range(1, 5)]
Rm = Ring('mod')
ch = FastChain(['A', 'F'], thN, hN, cvec, Rm)
vac = Rm.zeros(16); vac[3] = 1
rows = []
for a in range(4):
    b = ch.bra_t(vac, a, Rm.red(ch.th[0] - ch.h))
    for k in range(4):
        rows.append(b); b = ch.bra_t(b, 1, ch.th[1])
Wn = W.subs({th1: sp.Rational(2, 7), th2: sp.Rational(-3, 5), h: sp.Rational(1, 2)})
ok = all(int(rows[r][c]) == Rm.c(Fraction(int(sp.numer(Wn[r, c])), int(sp.denom(Wn[r, c])))) for r in range(16) for c in range(16))
print('fast mod-p code reproduces symbolic 2-site SoV matrix:', ok)
