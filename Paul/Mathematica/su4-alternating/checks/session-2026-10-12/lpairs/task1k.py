# Task 1k: (i) <W_a^L| (X3X1)^(1) = 24 h^7 y^2 l1 l3 z4^2 <W_{a z1/z4} (x) W_a^{L-1}|
#          (ii) pretwisted pair-p word: <W_{ar}^{p-1} (x) W_a^{L-p+1}| K^(p)_kappa = (1-k z1)^2 <W_{ar}^{p-1} (x) W_a' (x) W_a^{L-p}|
#          (iii) chain word (X3X1)^(1)...(X3X1)^(p-1) K^(p) on <W_a^L|
#          (iv) L = 4: pair-1 word exact; l = 2 multilinear span test
import numpy as np, sys, time
from fractions import Fraction as Fr
import eng, lp, modlin
from eng import P
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
r = z[0] / z[3]
for Lp in (3, 4):
    t0 = time.time()
    thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9), Fr(5, 6)][:Lp]
    ch = lp.chain_abjm(Lp, thF, h, z, 'p')
    a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
    for p in range(1, Lp + 1):
        uF = ch.th[2 * p - 1]
        l1, l3 = eng.lam(ch, 1, uF, p), eng.lam(ch, 3, uF + h, p)
        pre = lp.Wstate(ch, [a * r] * (p - 1) + [a] * (Lp - p + 1))
        x31 = ch.apply_t(1, uF, ch.apply_t(3, uF + h, pre, reg=True))
        tgt = eng.comb(ch, [24 * h**7 * y**2 * l1 * l3 * z[3]**2], [lp.Wstate(ch, [a * r] * p + [a] * (Lp - p))])
        res = [lp.eq(ch, x31, tgt)]
        for kap in (Fr(3, 10), 1 / z[3]):
            ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
            out = lp.applyK(ch, pre, p, kap, y)
            tgt = eng.comb(ch, [(1 - kap * z[0])**2], [lp.Wstate(ch, [a * r] * (p - 1) + [ap] + [a] * (Lp - p))])
            res.append(lp.eq(ch, out, tgt))
        print('L=%d p=%d: (X3X1)^(p) on pretwisted ref -> next pair twisted: %s ; K^(p) on pretwisted ref rotates pair p (kappa=3/10, 1/z4): %s'
              % (Lp, p, res[0], res[1:]))
    # (iii) chain word from the untwisted reference
    kap = Fr(3, 10); ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
    for p in range(2, Lp + 1):
        v = ref; c = Fr(1)
        for q in range(1, p):
            uF = ch.th[2 * q - 1]
            v = ch.apply_t(1, uF, ch.apply_t(3, uF + h, v, reg=True))
            c *= 24 * h**7 * y**2 * eng.lam(ch, 1, uF, q) * eng.lam(ch, 3, uF + h, q) * z[3]**2
        v = lp.applyK(ch, v, p, kap, y)
        tgt = eng.comb(ch, [c * (1 - kap * z[0])**2], [lp.Wstate(ch, [a * r] * (p - 1) + [ap] + [a] * (Lp - p))])
        print('   chain word (X3X1)^(1..%d) K^(%d) on <W_a^L| = c (1-kz1)^2 <W_ar^%d W_a\' W_a^%d| : %s' % (p - 1, p, p - 1, Lp - p, lp.eq(ch, v, tgt)))
    print('   (%.0f s)' % (time.time() - t0))
# (iv) L = 4, l = 2: multilinear span
Lp = 4; thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9), Fr(5, 6)]
ch = lp.chain_abjm(Lp, thF, h, z, 'p'); a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
t0 = time.time()
vecs = {(): ref}
for p in (1, 2):
    new = {}
    for key, v in vecs.items():
        new[key + (0,)] = v
        for i, x in enumerate(lp.Xs(ch, v, p)): new[key + (i + 1,)] = x
    vecs = new
mons = [vecs[k] for k in sorted(vecs)]
for kap in (Fr(3, 10), 1 / z[3]):
    ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
    tgt = lp.Wstate(ch, [ap, ap, a, a])
    rk, cons, sol, null = modlin.solve(mons, tgt)
    print('L=4 l=2 kappa=%s: 25 multilinear monomials, rank %d, target <W_a\'^2 W_a^2| in span: %s  (%.0f s)' % (kap, rk, cons, time.time() - t0))
