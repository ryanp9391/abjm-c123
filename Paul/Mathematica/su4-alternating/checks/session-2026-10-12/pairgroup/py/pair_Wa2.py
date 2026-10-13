# 1-pair chain at ABJM point, diag twist, symbolic: all localising words on W_a; solve for W_{a'}
import sympy as sp, numpy as np, itertools, pickle
from pgbase import *
h = sp.Symbol('h'); thF = sp.Integer(0); thA = 2*h
z = sp.symbols('z1:5'); a, kap, u = sp.symbols('a kappa u')
G = diagm(list(z))
ch = Chain(['A','F'], [thA, thF], h, G)
W = prod_state([[1,0,0,a],[-a,0,0,1]])
def T(k, m, psi):
    """t_k(thF + m h), regularised (vanishing scalar prefactors dropped, symbolic-safe)"""
    if k in (1,2) and not ((k == 2 and m in (0, 4))):
        return ch.apply_t(k, thF + m*h, psi)
    out = ch.apply_t(k, u, psi)
    pt = thF + m*h
    def red(e):
        e = sp.factor(sp.expand(e))
        if e == 0: return e
        while sp.simplify(e.subs(u, pt)) == 0:
            e = sp.cancel(e/(u-pt))
        return e.subs(u, pt)
    return np.vectorize(red, otypes=[object])(out)
words = [[], [(1,0)], [(2,1)], [(3,1)], [(1,0),(1,1)], [(1,0),(3,1)], [(1,0),(3,3)], [(1,2),(3,1)], [(3,1),(3,2)],
         [(1,0),(1,0),(3,3)], [(1,0),(1,2),(3,1)], [(1,0),(2,1),(3,3)], [(1,0),(2,2),(3,1)], [(1,0),(3,1),(3,3)], [(1,2),(2,1),(3,1)], [(1,2),(3,1),(3,1)]]
res = {}
for w in words:
    psi = W
    for (k, m) in w: psi = T(k, m, psi)
    psi = np.vectorize(lambda e: sp.factor(sp.expand(e)), otypes=[object])(psi)
    res[str(w)] = psi
    print(w, {(i+1,j+1): psi[i,j] for i in range(4) for j in range(4) if psi[i,j] != 0})
pickle.dump(res, open('pair_Wa_words.pkl','wb'))
