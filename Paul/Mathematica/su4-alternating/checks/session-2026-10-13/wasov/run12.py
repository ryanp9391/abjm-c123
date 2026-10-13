"""F-only bra coefficients at a second (h, theta) point: P_ab(rho2, rho4) for e1 (x) e_a (x) e1 (x) e_b, a,b in {2,3,4}; test the general-theta formulas."""
from wa_sov_core import *
import numpy as np
h = Fr(2, 3); th = (Fr(3, 11), Fr(-5, 7), Fr(9, 13), Fr(-13, 9))
ch = Chain(h=h, th=th)
W = sov_words(ch)
# inverse blockwise (only charges 0..4 needed for words (0,k2,0,k4), k<=2 -> charge <= 4)
Y = {}
for d, xs in charge_blocks().items():
    if d > 4: continue
    cs = [j for j, c in enumerate(COMPS) if comp_charge(c) == d]
    B = [[W[i][j] for j in cs] for i in xs]
    Bi = minv(B)
    for a, j in enumerate(cs):
        for b, i in enumerate(xs):
            Y[(j, i)] = Bi[a][b]
idxw = {w: i for i, w in enumerate(WORDS)}
def ci(l): return sum((c) * 4 ** (3 - p) for p, c in enumerate(l))
th12, th23, th24, th14, th34 = th[0] - th[1], th[1] - th[2], th[1] - th[3], th[0] - th[3], th[2] - th[3]
sc2 = h * th12 * th23 * (th24 - h); sc4 = h * th14 * th24 * th34
nu = th24
# formulas
pred = {
 (4, 4): {(0, 0): Fr(1)},
 (4, 3): {(0, 1): Fr(1), (1, 0): -h / nu},
 (3, 4): {(1, 0): Fr(1)},
 (3, 3): {(1, 1): nu / (nu + h)},
 (4, 2): {(0, 2): Fr(1), (1, 1): -(h / nu) * (1 + h / (nu + h)), (2, 0): -(h / nu) * (nu - h) / nu},
 (2, 4): {(2, 0): (nu - h) / nu, (1, 1): h / (nu + h)},
 (3, 2): {(1, 2): nu / (nu + h), (2, 1): -h / (nu + h)},
 (2, 3): {(2, 1): nu / (nu + h)},
 (2, 2): {(2, 2): (nu / (nu + h)) ** 2},
}
ok = True
for (a, b), poly in pred.items():
    for k2 in range(3):
        for k4 in range(3):
            c = Y.get((ci([0, a - 1, 0, b - 1]), idxw[(0, k2, 0, k4)]), Fr(0)) * sc2 ** k2 * sc4 ** k4
            want = poly.get((k2, k4), Fr(0))
            if c != want:
                ok = False; print('mismatch', (a, b), (k2, k4), c, want)
print('general-theta F-only formulas reproduced at (h,theta) = (2/3, (3/11,-5/7,9/13,-13/9)):', ok)
