# exact (Fraction) 4-site check: span of localising words on W_a^2 at the ABJM point vs W_{a'} (x) W_a
import numpy as np, sympy as sp
from pgbase import *
h = Fr(1,2); thF = Fr(-3,5); thA = thF + 2*h
th = [thA, thF, Fr(-3,8), Fr(-11,8)]
z = [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
ch = Chain(['A','F','A','F'], th, h, diagm(z))
a = Fr(1,3)
Wp = lambda a: [[Fr(1),Fr(0),Fr(0),a],[-a,Fr(0),Fr(0),Fr(1)]]
ref = prod_state(Wp(a) + Wp(a))
rest = prod_state(Wp(a)).reshape(-1); rr = rest.dot(rest)
words = [[], [(1,0)], [(2,1)], [(3,1)], [(1,0),(1,1)], [(1,0),(3,1)], [(1,0),(3,3)], [(1,2),(3,1)], [(3,1),(3,2)], [(1,0),(1,0),(3,3)]]
cols = []
for w in words:
    psi = ref
    for (k,m) in w: psi = ch.apply_t(k, thF + m*h, psi, reg=True)
    M = psi.reshape(16, 16); ps = M.dot(rest)/rr
    assert all(x == 0 for x in (M - np.multiply.outer(ps, rest)).reshape(-1)), w
    cols.append(ps)
B = sp.Matrix(16, len(cols), lambda i, j: cols[j][i])
print('exact rank of localised span on W_a^2:', B.rank())
kap = Fr(3,10)
ap = a*(1-kap*z[3])/(1-kap*z[0])
tgt = sp.Matrix(prod_state(Wp(ap)).reshape(-1).tolist())
aug = B.row_join(tgt)
print('rank with target W_{a\'}:', aug.rank())
