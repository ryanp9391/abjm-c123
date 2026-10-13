# Task 1b: naive l-pair products prod_{p<=l} K^(p) on <W_a^L| (generic pair-preserving theta, ABJM point), exact mod p.
# Also: support analysis of <W_a^L| K^(p) for single p (which pairs are changed).
import numpy as np, itertools
from fractions import Fraction as Fr
import eng, lp
from eng import P
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]

def support_pairs(ch, vec, alist):
    """smallest set S of pairs such that vec = (something on pairs S) (x) (W_{a_n} on pairs not in S); test by projection mod p"""
    Lp = ch.Lp
    for r in range(0, Lp + 1):
        for S in itertools.combinations(range(Lp), r):
            # contract vec with the dual of W_{a_n} on pairs outside S: check vec is in the image
            # build projector test: vec reshaped pairs-> (16,)*Lp ; outside S must be proportional to W_{a_n}
            V = vec.reshape((16,) * Lp)
            ok = True
            for n in range(Lp):
                if n in S: continue
                w = eng.prod_state(eng.W(alist[n]), 'p').reshape(16)
                # vec must be (rest) (x) w on pair n: check V_{..i..} * w_j == V_{..j..} * w_i for all i,j via a pivot
                piv = int(np.nonzero(w)[0][0])
                Vn = np.moveaxis(V, n, -1)
                lhs = (Vn[..., :] * w[piv]) % P
                rhs = (Vn[..., piv:piv + 1] * w[None, :]) % P
                if not np.array_equal(lhs, rhs): ok = False; break
            if ok: return S
    return None

for Lp in (2, 3):
    thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)][:Lp]
    ch = lp.chain_abjm(Lp, thF, h, z, 'p')
    a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
    kap = Fr(3, 10); ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
    for p in range(1, Lp + 1):
        out = lp.applyK(ch, ref, p, kap, y)
        print('L=%d  <W_a^L|K^(%d>: changed pairs (0-based) = %s' % (Lp, p, support_pairs(ch, out, [a] * Lp)))
        for i, x in enumerate(lp.Xs(ch, ref, p)):
            print('      word %s of pair %d: changed pairs %s' % (['X1', 'X2', 'X3', 'X3X1'][i], p, support_pairs(ch, x, [a] * Lp)))
    for l in range(1, Lp + 1):
        for kap in (Fr(3, 10), Fr(-7, 5), 1 / z[3]):
            ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
            v = ref
            for p in range(1, l + 1): v = lp.applyK(ch, v, p, kap, y)
            tgt = eng.comb(ch, [(1 - kap * z[0])**(2 * l)], [lp.Wstate(ch, [ap] * l + [a] * (Lp - l))])
            print('L=%d l=%d kappa=%s: naive product = (1-kz1)^(2l) <W_a\'^l W_a^(L-l)| : %s ; changed pairs of result: %s'
                  % (Lp, l, kap, lp.eq(ch, v, tgt), support_pairs(ch, v, [a] * Lp)))
