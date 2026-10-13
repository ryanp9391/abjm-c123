# zero-weight blocks (E11,E22,E33,E44) of pair-chain operators at the ABJM point, diagonal twist, symbolic z, h
import numpy as np, sympy as sp
from pgbase import *
h = sp.Symbol('h'); thF = sp.Integer(0); thA = 2*h
z = sp.symbols('z1:5')
ch = Chain(['A','F'], [thA, thF], h, diagm(list(z)))
zw = [(0,0),(1,1),(2,2),(3,3)]
def block(k, u0):
    B = sp.zeros(4,4)
    for r,(i,i_) in enumerate(zw):
        e = np.zeros((4,4), dtype=object); e[:] = 0; e[i,i] = 1
        out = ch.apply_t(k, u0, e, reg=True)
        for c,(j,j_) in enumerate(zw): B[r,c] = sp.factor(sp.expand(out[j,j]))
    return B
X1 = block(1, thF); X2 = block(2, thF+h); X3 = block(3, thF+h); Y = block(1, thF+h); T1A = block(1, thA)
for n, B in (('X1 t1(F)', X1), ('X2 t2(F+h)', X2), ('X3 t3^(F+h)', X3), ('t1(F+h)', Y), ('t1(A)', T1A)):
    print(n); sp.pprint(B, wrap_line=False)
# algebra spanned by localising words on the zero-weight block
mats = {'1': sp.eye(4), 'X1': X1, 'X2': X2, 'X3': X3, 'X3X1': X3*X1}
M = sp.Matrix([list(m.reshape(16,1)) for m in mats.values()])
print('rank of span{1,X1,X2,X3,X3X1} on zero-weight block:', M.rank())
M2 = sp.Matrix([list(m.reshape(16,1)) for m in list(mats.values()) + [Y]])
print('... adding t1(F+h):', M2.rank())
