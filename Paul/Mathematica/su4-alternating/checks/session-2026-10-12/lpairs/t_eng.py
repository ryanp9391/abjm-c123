# validate eng.py against pgbase.py (exact Fractions) at 2 pairs, random product bra, k = 1,2,3 (+hat)
import numpy as np, time
from fractions import Fraction as Fr
import pgbase, eng
h = Fr(1, 2); th = [Fr(2, 7), Fr(-3, 5), Fr(8, 11), Fr(-11, 8)]; z = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
rng = np.random.default_rng(1)
vecs = [[Fr(int(x), 7) for x in rng.integers(-5, 6, 4)] for _ in range(4)]
ref = pgbase.prod_state(vecs)
chP = pgbase.Chain(['A', 'F'] * 2, th, h, pgbase.diagm(z))
chp = eng.Chain(2, th, h, z, 'p'); chc = eng.Chain(2, [float(x) for x in th], 0.5, [float(x) for x in z], 'c')
refp = eng.prod_state(vecs, 'p'); refc = eng.prod_state([[float(x) for x in v] for v in vecs], 'c')
for k, u, reg in [(1, Fr(3, 10), False), (2, Fr(-7, 9), False), (3, Fr(5, 13), False), (3, th[1] + h, True), (2, th[1] + h, False)]:
    ex = chP.apply_t(k, u, ref, reg=reg)
    exp_ = np.vectorize(eng.modp, otypes=[np.int64])(ex)
    outp = chp.apply_t(k, u, refp, reg=reg)
    outc = chc.apply_t(k, float(u), refc, reg=reg)
    exf = ex.astype(float)
    print(k, u, reg, 'mod-p equal:', np.array_equal(exp_, outp), ' complex rel err: %.1e' % (np.abs(outc - exf).max() / np.abs(exf).max()))
# commutativity check t1(u) t2(v) on chain, mod p
a = chp.apply_t(2, Fr(1, 3), chp.apply_t(1, Fr(2, 9), refp)); b = chp.apply_t(1, Fr(2, 9), chp.apply_t(2, Fr(1, 3), refp))
print('[t1,t2] = 0 mod p:', np.array_equal(a, b))
# timing at L = 4 (65536)
ch4 = eng.Chain(4, [Fr(i, 7) for i in range(8)], h, z, 'p')
r4 = eng.prod_state(vecs * 2, 'p')
for k in (1, 2, 3):
    t0 = time.time(); ch4.apply_t(k, Fr(1, 3), r4); print('L=4 mod-p t%d: %.2f s' % (k, time.time() - t0))
ch4c = eng.Chain(4, [i / 7 for i in range(8)], 0.5, [1.5, 1.25, 2 / 3, 4 / 7], 'c')
r4c = eng.prod_state([[float(x) for x in v] for v in vecs * 2], 'c')
for k in (1, 2, 3):
    t0 = time.time(); ch4c.apply_t(k, 1 / 3, r4c); print('L=4 complex t%d: %.2f s' % (k, time.time() - t0))
