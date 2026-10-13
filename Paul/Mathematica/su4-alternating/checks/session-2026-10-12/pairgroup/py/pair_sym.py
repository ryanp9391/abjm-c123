# symbolic: <e1 e4| X_k on the 1-pair chain at the ABJM point, generic twist G
import sympy as sp, numpy as np
from pgbase import *
h = sp.Symbol('h'); thF = sp.Integer(0); thA = 2*h
Gs = sp.Matrix(4, 4, lambda i, j: sp.Symbol('g%d%d' % (i+1, j+1)))
G = np.array(Gs.tolist(), dtype=object)
ch = Chain(['A','F'], [thA, thF], h, G)
vac = np.zeros((4,4), dtype=object); vac[:] = 0; vac[0,3] = 1
u = sp.Symbol('u')
def show(name, out):
    M = sp.Matrix(4,4, lambda i,j: sp.factor(sp.expand(out[i,j])))
    print(name); sp.pprint(M, wrap_line=False); print()
show('X1 = t1(thF)', ch.apply_t(1, thF, vac))
show('X2 = t2(thF+h)', ch.apply_t(2, thF + h, vac))
t3u = ch.apply_t(3, u, vac)
X3 = np.vectorize(lambda e: sp.factor(sp.cancel(sp.expand(e)/(u - thF - h))).subs(u, thF + h), otypes=[object])(t3u)
show('X3 = t3hat(thF+h)', X3)
show('t4(u)', ch.apply_t(4, u, vac))
