# 1-pair chain at the ABJM point, diagonal twist: eigen-decomposition; eigenvalues of localising words on the
# zero-weight eigenstates; which pair-chain operators separate them
import numpy as np, sympy as sp
from pgbase import *
h = Fr(1,2); thF = Fr(0); thA = 2*h
z = [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
ch = Chain(['A','F'], [thA, thF], h, diagm(z))
def opm(k, u0, reg=True):
    M = sp.zeros(16,16)
    for i in range(16):
        e = np.zeros(16, dtype=object); e[:] = Fr(0); e[i] = Fr(1)
        out = ch.apply_t(k, u0, e.reshape(4,4), reg=reg).reshape(-1)
        for j in range(16): M[i,j] = sp.Rational(out[j].numerator, out[j].denominator)
    return M
X1 = opm(1, thF); X2 = opm(2, thF+h); X3 = opm(3, thF+h); Y1 = opm(1, thF+h); T1A = opm(1, thA); T33 = opm(3, thF+3*h); T32 = opm(3, thF+2*h)
gen = opm(1, Fr(37,100)) + 3*opm(3, Fr(-29,100))  # generic element to diagonalise
# zero-weight subspace indices: E_ii -> index 4*i+i
zw = [0, 5, 10, 15]
def restrict(M): return M.extract(zw, zw)
print('zero-weight block of a generic t (bra action, rows = input):')
Gz = restrict(gen)
# bra eigenvectors: psi M = lam psi  -> left eigenvectors = eigenvectors of M^T
ev = (Gz.T).eigenvects()
print('number of zero-weight eigenvalues:', len(ev))
vecs = [v[2][0].T for v in ev]
for name, M in (('X1=t1(F)', X1), ('X2=t2(F+h)', X2), ('X3=t3^(F+h)', X3), ('t1(F+h)', Y1), ('t1(A)', T1A), ('t3^(F+3h)', T33), ('t3(F+2h)=t3(A)', T32)):
    Mz = restrict(M)
    vals = []
    for v in vecs:
        w = v*Mz
        # ratio
        idx = [i for i in range(4) if v[i] != 0][0]
        lam = sp.nsimplify(w[idx]/v[idx])
        vals.append(sp.N(lam, 8))
        assert sp.simplify(w - lam*v) == sp.zeros(1,4), name
    print(name, vals)
