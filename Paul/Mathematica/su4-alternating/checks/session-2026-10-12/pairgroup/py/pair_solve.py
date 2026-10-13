# Pair level (1-pair chain), ABJM point, diagonal twist with z2 = z3 = y: solve
#   sum_w c_w <W_a| w = N(kappa) <W_{a'}|,  a' = a (1-kappa z4)/(1-kappa z1),  w in {1, X1, X2, X3, X3X1}
import sympy as sp, numpy as np, pickle
from pgbase import *
h = sp.Symbol('h', positive=True); thF = sp.Integer(0); thA = 2*h
z1, y, z4, a, kap = sp.symbols('z1 y z4 a kappa')
zz = [z1, y, y, z4]
ch = Chain(['A','F'], [thA, thF], h, diagm(zz))
W = prod_state([[1,0,0,a],[-a,0,0,1]])
X = {}
X['1'] = W
X['X1'] = ch.apply_t(1, thF, W)
X['X2'] = ch.apply_t(2, thF+h, W)
X['X3'] = ch.apply_t(3, thF+h, W, reg=True)
X['X3X1'] = ch.apply_t(1, thF, X['X3'])
names = list(X)
cs = sp.symbols('c0:5')
N = (1 - kap*z1)**2          # normalisation making the target polynomial in kappa
ap = a*(1-kap*z4)/(1-kap*z1)
tgt = N*prod_state([[1,0,0,ap],[-ap,0,0,1]])
eqs = []
for i in range(4):
    for j in range(4):
        e = sum(c*X[n][i,j] for c, n in zip(cs, names)) - tgt[i,j]
        e = sp.together(sp.expand(e))
        if e != 0: eqs.append(sp.numer(e))
sol = sp.solve(eqs, cs, dict=True)
print(len(sol))
for s in sol:
    for c in cs: print(c, '=', sp.factor(s[c]))
pickle.dump({str(k): v for k, v in sol[0].items()}, open('pair_solution.pkl','wb'))
