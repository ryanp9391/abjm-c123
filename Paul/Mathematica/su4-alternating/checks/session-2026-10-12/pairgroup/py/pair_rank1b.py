# rank-1 elements of span{E14, X1, X2, X3, X3X1} on the 1-pair chain at the ABJM point, generic rational G
import sympy as sp, numpy as np, random
from pgbase import *
random.seed(5)
h = sp.Rational(1,2); thF = sp.Integer(0); thA = 2*h
Gs = sp.Matrix(4,4, lambda i,j: sp.Rational(random.randint(-9,9), random.randint(1,5)))
G = np.array(Gs.tolist(), dtype=object)
ch = Chain(['A','F'], [thA, thF], h, G)
vac = np.zeros((4,4), dtype=object); vac[:] = 0; vac[0,3] = 1
u = sp.Symbol('u')
def t3hat(psi):
    t3u = ch.apply_t(3, u, psi)
    return np.vectorize(lambda e: sp.cancel(sp.expand(e)/(u-thF-h)).subs(u, thF+h), otypes=[object])(t3u)
X1 = ch.apply_t(1, thF, vac); X2 = ch.apply_t(2, thF+h, vac); X3 = t3hat(vac); X31 = ch.apply_t(1, thF, X3)
c = sp.symbols('c1:5')
M = sp.Matrix(vac.tolist()) + c[0]*sp.Matrix(X1.tolist()) + c[1]*sp.Matrix(X2.tolist()) + c[2]*sp.Matrix(X3.tolist()) + c[3]*sp.Matrix(X31.tolist())
eqs = set()
for i in range(4):
    for j in range(i+1,4):
        for k in range(4):
            for l in range(k+1,4):
                e = sp.expand(M[i,k]*M[j,l]-M[i,l]*M[j,k])
                if e != 0: eqs.add(e)
sol = sp.solve(list(eqs), list(c), dict=True)
print(sol)
