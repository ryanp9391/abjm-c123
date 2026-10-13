# Task 1a: reproduce the pair-1 word; single pair-p words K^(p) on <W_a^L|; coordinator's identity t2(u0) = t1(u0) t1(thF)
import numpy as np, sys, time
from fractions import Fraction as Fr
import eng, lp
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
for geom in ('generic', 'periodic'):
    for Lp in (2, 3):
        thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)][:Lp] if geom == 'generic' else [Fr(-3, 5)] * Lp
        ch = lp.chain_abjm(Lp, thF, h, z, 'p')
        a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
        # identity checks on a random product bra
        rng = np.random.default_rng(5)
        rb = eng.prod_state([[Fr(int(v), 3) for v in rng.integers(-4, 5, 4)] for _ in range(2 * Lp)], 'p')
        u0 = thF[0] + h
        id1 = lp.eq(ch, ch.apply_t(2, u0, rb), ch.apply_t(1, u0, ch.apply_t(1, thF[0], rb)))
        id2 = not ch.apply_t(3, u0, rb).any()
        print('%s L=%d: t2(u0) = t1(u0) t1(thF) on random bra: %s;  t3(u0) (full) = 0: %s' % (geom, Lp, id1, id2))
        for kap in (Fr(3, 10), Fr(-7, 5), 1 / z[3]):
            ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
            res = []
            for p in range(1, Lp + 1):
                out = lp.applyK(ch, ref, p, kap, y)
                al = [a] * Lp; al[p - 1] = ap
                tgt = eng.comb(ch, [(1 - kap * z[0])**2], [lp.Wstate(ch, al)])
                res.append(lp.eq(ch, out, tgt))
            print('   kappa=%s: <W_a^L|K^(p) = (1-kz1)^2 <..W_a\'(pair p)..|  for p=1..%d: %s' % (kap, Lp, res))
