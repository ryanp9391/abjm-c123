# Task 1g: right-end localisation scan at L = 3 (ABJM point, W_a^L, z2 = z3): t_k(u) for u = th_{A,L}+mh, th_{F,L}+mh
from fractions import Fraction as Fr
import eng, lp
from supp import support_pairs
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)]
ch = lp.chain_abjm(3, thF, h, z, 'p')
a = Fr(1, 3); ref = lp.Wstate(ch, [a] * 3)
for site, base in (('A3', ch.th[4]), ('F3', ch.th[5]), ('A2', ch.th[2]), ('F2', ch.th[3])):
    for k in (1, 2, 3):
        for m in range(-4, 5):
            u = base + m * h
            out = ch.apply_t(k, u, ref, reg=True)
            if not out.any(): print('  t%d(th_%s%+dh): ZERO' % (k, site, m)); continue
            S = support_pairs(ch, out, [a] * 3)
            if S is not None and len(S) < 3:
                print('  t%d(th_%s%+dh): changed pairs %s' % (k, site, m, S))
