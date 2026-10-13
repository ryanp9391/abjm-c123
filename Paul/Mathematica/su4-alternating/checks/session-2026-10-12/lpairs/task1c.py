# Task 1c: is the target <W_a'^l (x) W_a^(L-l)| in the span of the 5^l monomials prod_{p<=l} w_p, w_p in {1,X1,X2,X3,X3X1}^(p) ?
import numpy as np, itertools, time
from fractions import Fraction as Fr
import eng, lp, modlin
from eng import P
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
names = ['1', 'X1', 'X2', 'X3', 'X3X1']

def monomials(ch, ref, l):
    vecs = {(): ref}
    for p in range(1, l + 1):
        new = {}
        for key, v in vecs.items():
            new[key + (0,)] = v
            for i, x in enumerate(lp.Xs(ch, v, p)): new[key + (i + 1,)] = x
        vecs = new
    keys = sorted(vecs)
    return keys, [vecs[k] for k in keys]

for Lp in (2, 3):
    thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)][:Lp]
    ch = lp.chain_abjm(Lp, thF, h, z, 'p')
    a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
    for l in range(2, Lp + 1):
        t0 = time.time()
        keys, mons = monomials(ch, ref, l)
        for kap in (Fr(3, 10), 1 / z[3]):
            ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
            tgt = eng.comb(ch, [(1 - kap * z[0])**(2 * l)], [lp.Wstate(ch, [ap] * l + [a] * (Lp - l))])
            rk, cons, sol, null = modlin.solve(mons, tgt)
            print('L=%d l=%d kappa=%s: %d monomials, rank %d, target in span: %s, nullity %d  (%.1f s)'
                  % (Lp, l, kap, len(mons), rk, cons, len(null), time.time() - t0))
