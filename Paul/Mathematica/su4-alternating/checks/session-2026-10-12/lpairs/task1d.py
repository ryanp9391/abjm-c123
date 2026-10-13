# Task 1d: Krylov closure of <W_a^L| under the pair-point operators of pairs 1..l (t1(thF_p), t2(thF_p+h), t3^(thF_p+h)),
# L = 3; is the target <W_a'^l W_a^(L-l)| in the closure?  Also closure under pair-1 ops only (should be 5).
import numpy as np, time
from fractions import Fraction as Fr
import eng, lp
from modspan import Span
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
Lp = 3; thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)]
ch = lp.chain_abjm(Lp, thF, h, z, 'p')
a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
def ops(pairs):
    out = []
    for p in pairs:
        uF = ch.th[2 * p - 1]
        out += [lambda v, u=uF: ch.apply_t(1, u, v), lambda v, u=uF: ch.apply_t(2, u + h, v), lambda v, u=uF: ch.apply_t(3, u + h, v, reg=True)]
    return out
for pairs in ([1], [1, 2], [1, 2, 3]):
    t0 = time.time()
    S = Span(4 ** (2 * Lp)); S.add(ref); frontier = [ref]; G = ops(pairs)
    while frontier:
        new = []
        for v in frontier:
            for g in G:
                w = g(v)
                if S.add(w): new.append(w)
        frontier = new
    kap = Fr(3, 10); ap = a * (1 - kap * z[3]) / (1 - kap * z[0]); l = len(pairs)
    tgt = lp.Wstate(ch, [ap] * l + [a] * (Lp - l))
    tgt_end = lp.Wstate(ch, [a] * (Lp - l) + [ap] * l)
    print('pairs %s: closure dim %d, target (first l pairs rotated) in closure: %s, (last l pairs rotated): %s  (%.1f s)'
          % (pairs, S.dim(), S.contains(tgt), S.contains(tgt_end), time.time() - t0))
