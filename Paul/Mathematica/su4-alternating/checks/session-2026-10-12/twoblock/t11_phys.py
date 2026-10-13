"""Physical HW class: B = (W_{-1} x vac) rho(n D), with pair 2 = companion SoV vacuum.
Diag frame: pair1 = (lam - g e4) x (1 + (g/lam1) e1), pair2 = lam x 1 ;  companion: pair2 = e1 x e4."""
import numpy as np
from core import *
from exact import *
from fractions import Fraction as Fr
import sympy
R8 = res8(); ps = list(R8.keys())
z = z0
Vs = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in Vand(z)])
Vi = Vs.inv()
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
print('lambda_j =', lam, ' sum =', sum(lam))
def tofr(x): return Fr(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1]))
VT = [[tofr(Vs[j, i]) for j in range(4)] for i in range(4)]
Vinv = [[tofr(Vi[i, j]) for j in range(4)] for i in range(4)]
def vm(v, A): return [sum(v[i] * A[i][j] for i in range(4)) for j in range(4)]
def A2c(v): return vm(v, VT)      # A site diag -> companion
def F2c(u): return vm(u, Vinv)    # F site diag -> companion
one = [Fr(1)] * 4
print('companion image of lam (A):', A2c(lam), '  of (1,1,1,1) (F):', F2c(one))
def coeffs_exact(vs):
    """vs: 4 companion-frame site vectors (Fractions). returns dict label -> Fraction"""
    out = {}
    per_p = {}
    for p in ps:
        R = Ring('mod', p)
        b = kronvecs(R, [R.arr(v) for v in vs])
        per_p[p] = R.mm(b[None, :], R8[p][1])[0]
    for x in range(256):
        if per_p[ps[0]][x] == 0 and per_p[ps[1]][x] == 0: continue
        out[wordsIdx[x]] = reconstruct({p: int(per_p[p][x]) for p in ps})
    return out
e1, e4 = [1, 0, 0, 0], [0, 0, 0, 1]
# gamma-expansion: pair1(diag) = lam x 1  - g [ e4 x 1 - lam x e1/lam1 ] - g^2 e4 x e1/lam1
terms = {
  'g^0': [A2c(lam), F2c(one)],
  'g^1 (a)': [A2c(e4), F2c(one)],          # coefficient -1
  'g^1 (b)': [A2c(lam), F2c([1 / lam[0], 0, 0, 0])],  # coefficient +1
  'g^2': [A2c(e4), F2c([1 / lam[0], 0, 0, 0])],      # coefficient -1
}
pair2 = [A2c(lam), F2c(one)]
C = {}
for k, pr in terms.items():
    C[k] = coeffs_exact(pr + pair2)
    print(k, 'nnz', len(C[k]), ' failures', sum(v is None for v in C[k].values()))
import pickle; pickle.dump((C, lam), open('phys_coeffs.pkl', 'wb'))
for k in C:
    print('---', k)
    for x, v in sorted(C[k].items()):
        print(x, v)
