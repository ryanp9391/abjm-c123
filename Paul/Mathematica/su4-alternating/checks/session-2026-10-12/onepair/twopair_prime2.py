# rerun the two-pair factorisation tests with a second prime (robustness of mod-p ranks)
import numpy as np, random
from fractions import Fraction
from itertools import combinations
import fastchain as fc
fc.P = 1000000007
P = fc.P
hN = Fraction(1, 2); thN = [Fraction(2, 7), Fraction(-3, 5), Fraction(8, 11), Fraction(-11, 8)]
zA = [Fraction(3, 2), Fraction(5, 4), Fraction(2, 3), Fraction(4, 7)]
cv = [sum(np.prod(s) for s in combinations(zA, k)) for k in range(1, 5)]
R = fc.Ring('mod')
ch = fc.FastChain(['A', 'F', 'A', 'F'], thN, hN, cv, R)
vac = R.zeros(256); vac[((0 * 4 + 3) * 4 + 0) * 4 + 3] = 1
rows = {}
for a1 in range(4):
    b1 = ch.bra_t(vac, a1, R.red(ch.th[0] - ch.h))
    for k2 in range(4):
        for a3 in range(4):
            b3 = ch.bra_t(b1, a3, R.red(ch.th[2] - ch.h))
            for k4 in range(4):
                rows[(a1, k2, a3, k4)] = b3; b3 = ch.bra_t(b3, 1, ch.th[3])
        b1 = ch.bra_t(b1, 1, ch.th[1])
labs = sorted(rows); W = np.array([rows[l] for l in labs], dtype=object)
print('rank W (p2):', fc.rank_mod(W))
Winv = fc.inv_mod(W)
def C_of(v1, u1, v2, u2):
    x = np.kron(np.kron(np.kron(np.array(v1, dtype=object), np.array(u1, dtype=object)), np.array(v2, dtype=object)), np.array(u2, dtype=object)) % P
    c = x.dot(Winv) % P
    C = np.zeros((16, 16), dtype=object)
    for n, (a1, k2, a3, k4) in enumerate(labs): C[4 * a1 + k2, 4 * a3 + k4] = c[n]
    return C
random.seed(99)
rnd = lambda: [random.randrange(P) for _ in range(4)]
e = lambda i: [1 if k == i else 0 for k in range(4)]
for name, vs in [('generic', (rnd(), rnd(), rnd(), rnd())), ('pair1 vac', (e(0), e(3), rnd(), rnd())), ('pair2 vac', (rnd(), rnd(), e(0), e(3)))]:
    print(name, 'rank', fc.rank_mod(C_of(*vs)))
from collections import Counter
print('ranks of coefficient bilinear forms:', sorted(Counter(fc.rank_mod(Winv[:, n].reshape(16, 16)) for n in range(256)).items()))
