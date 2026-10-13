# find c's s.t. E14 + c1 X1 + c2 X2 + c3 X3 + c4 detG*1 is rank one (product pair state); generic rational G
import sympy as sp, numpy as np, random
from pgbase import *
random.seed(5)
h = sp.Rational(1,2); thF = sp.Integer(0); thA = 2*h
Gs = sp.Matrix(4,4, lambda i,j: sp.Rational(random.randint(-9,9), random.randint(1,5)))
G = np.array(Gs.tolist(), dtype=object)
ch = Chain(['A','F'], [thA, thF], h, G)
vac = np.zeros((4,4), dtype=object); vac[:] = 0; vac[0,3] = 1
u = sp.Symbol('u')
X1 = sp.Matrix(ch.apply_t(1, thF, vac).tolist())
X2 = sp.Matrix(ch.apply_t(2, thF+h, vac).tolist())
t3u = ch.apply_t(3, u, vac)
X3 = sp.Matrix(4,4, lambda i,j: sp.cancel(sp.expand(t3u[i,j])/(u-thF-h)).subs(u, thF+h))
dG = Gs.det()
c1,c2,c3,c4 = sp.symbols('c1:5')
E14 = sp.zeros(4,4); E14[0,3] = 1
M = E14 + c1*X1 + c2*X2 + c3*X3 + c4*dG*sp.eye(4)
eqs = []
for i in range(4):
    for j in range(i+1,4):
        for k in range(4):
            for l in range(k+1,4):
                eqs.append(sp.expand(M[i,k]*M[j,l]-M[i,l]*M[j,k]))
eqs = [e for e in eqs if e != 0]
print(len(eqs))
sol = sp.solve(eqs, [c1,c2,c3,c4], dict=True)
print(sol)
