# 1-pair chain at ABJM point, diag twist, symbolic: all localising words on W_a (reg = drop vanishing scalar factors)
import sympy as sp, numpy as np, itertools, pickle
from pgbase import *
h = sp.Symbol('h'); thF = sp.Integer(0); thA = 2*h
z = sp.symbols('z1:5'); a, kap = sp.symbols('a kappa')
G = diagm(list(z))
ch = Chain(['A','F'], [thA, thF], h, G)
W = prod_state([[1,0,0,a],[-a,0,0,1]])
words = [[], [(1,0)], [(2,1)], [(3,1)], [(1,0),(1,1)], [(1,0),(3,1)], [(1,0),(3,3)], [(1,2),(3,1)], [(3,1),(3,2)],
         [(1,0),(1,0),(3,3)], [(1,0),(1,2),(3,1)], [(1,0),(2,1),(3,3)], [(1,0),(2,2),(3,1)], [(1,0),(3,1),(3,3)], [(1,2),(2,1),(3,1)], [(1,2),(3,1),(3,1)]]
res = {}
for w in words:
    psi = W
    for (k, m) in w: psi = ch.apply_t(k, thF + m*h, psi, reg=True)
    psi = np.vectorize(lambda e: sp.factor(sp.expand(e)), otypes=[object])(psi)
    res[str(w)] = psi
    print(w, {(i+1,j+1): psi[i,j] for i in range(4) for j in range(4) if psi[i,j] != 0})
pickle.dump(res, open('pair_Wa_words.pkl','wb'))
