# Task 1e: block localisation of monomials in the pair-point operators of pairs 1,2:
#   <W_a^3| m  =?=  lambda_m * (<W_a^2| m on the 2-pair chain) (x) <W_a|     (exact mod p)
import numpy as np, itertools, time
from fractions import Fraction as Fr
import eng, lp
from eng import P, modp
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)]
ch3 = lp.chain_abjm(3, thF, h, z, 'p'); ch2 = lp.chain_abjm(2, thF[:2], h, z, 'p')
a = Fr(1, 3); r3 = lp.Wstate(ch3, [a] * 3); r2 = lp.Wstate(ch2, [a] * 2)
wa = eng.prod_state(eng.W(a), 'p').reshape(16)
def gen(ch, v, g):
    p, k = g; uF = ch.th[2 * p - 1]
    if k == 1: return ch.apply_t(1, uF, v)
    if k == 2: return ch.apply_t(2, uF + h, v)
    return ch.apply_t(3, uF + h, v, reg=True)
G = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3)]
names = {(1, 1): 'X1a', (1, 2): 'X2a', (1, 3): 'X3a', (2, 1): 'X1b', (2, 2): 'X2b', (2, 3): 'X3b'}
def blockloc(v3, v2):
    V = v3.reshape(256, 16) % P
    # v3 = (v2 * c) (x) wa ?
    W2 = v2.reshape(256) % P
    if not W2.any(): return 'pair-chain result 0'
    i = int(np.nonzero(W2)[0][0]); j = int(np.nonzero(wa)[0][0])
    c = V[i, j] * pow(int(W2[i] * wa[j] % P), P - 2, P) % P
    ok = np.array_equal(V, (c * (np.multiply.outer(W2, wa) % P)) % P)
    return ('YES lam=%s' % c) if ok else 'no'
for deg in (1, 2):
    for mono in itertools.combinations_with_replacement(G, deg):
        v3, v2 = r3, r2
        for g in mono: v3 = gen(ch3, v3, g); v2 = gen(ch2, v2, g)
        print(deg, '*'.join(names[g] for g in mono), blockloc(v3, v2))
