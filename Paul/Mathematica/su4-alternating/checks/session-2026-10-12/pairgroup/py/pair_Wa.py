# 1-pair chain at the ABJM point, diagonal twist: <W_a| X for the 5 localising pair words; solve for W_{a'}
import sympy as sp, numpy as np
from pgbase import *
h = sp.Symbol('h'); thF = sp.Integer(0); thA = 2*h
z = sp.symbols('z1:5'); a, kap, u = sp.symbols('a kappa u')
G = diagm(list(z))
ch = Chain(['A','F'], [thA, thF], h, G)
def Wpair(a):
    v = [1,0,0,a]; w = [-a,0,0,1]
    return prod_state([v, w])
def t3hat(psi):
    t3u = ch.apply_t(3, u, psi)
    return np.vectorize(lambda e: sp.cancel(sp.expand(e)/(u-thF-h)).subs(u, thF+h), otypes=[object])(t3u)
simp = np.vectorize(lambda e: sp.factor(sp.expand(e)), otypes=[object])
W = Wpair(a)
X = {'1': W,
     'X1=t1(F)': ch.apply_t(1, thF, W),
     'X2=t2(F+h)': ch.apply_t(2, thF+h, W),
     'X3=t3^(F+h)': t3hat(W)}
X['X3X1'] = ch.apply_t(1, thF, X['X3=t3^(F+h)'])
for k, v in X.items():
    print(k); sp.pprint(sp.Matrix(simp(v).tolist()), wrap_line=False); print()
