import numpy as np, pickle, time
from fractions import Fraction
from itertools import combinations
from fastchain import *
hN = Fraction(1, 2); thN = [Fraction(2, 7), Fraction(-3, 5), Fraction(8, 11), Fraction(-11, 8)]
def cv(z): return [sum(np.prod(s) for s in combinations(z, k)) for k in range(1, 5)]
zA = [Fraction(3, 2), Fraction(5, 4), Fraction(2, 3), Fraction(4, 7)]
zB = [Fraction(7, 3), Fraction(-2, 5), Fraction(5, 4), Fraction(3)]
def basis(z, ring):
    ch = FastChain(['A', 'F', 'A', 'F'], thN, hN, cv(z), ring)
    vac = ring.zeros(256); vac[((0 * 4 + 3) * 4 + 0) * 4 + 3] = 1
    rows = {}
    for a1 in range(4):
        b1 = ch.bra_t(vac, a1, ring.red(ch.th[0] - ch.h))
        for k2 in range(4):
            for a3 in range(4):
                b3 = ch.bra_t(b1, a3, ring.red(ch.th[2] - ch.h))
                for k4 in range(4):
                    rows[(a1, k2, a3, k4)] = b3
                    b3 = ch.bra_t(b3, 1, ch.th[3])
            b1 = ch.bra_t(b1, 1, ch.th[1])
    return rows
t0 = time.time()
Rm = Ring('mod')
rowsA = basis(zA, Rm); print('basis z built', time.time() - t0)
rowsB = basis(zB, Rm); print('basis z\' built', time.time() - t0)
labs = sorted(rowsA)
WA = np.array([rowsA[l] for l in labs], dtype=object); WB = np.array([rowsB[l] for l in labs], dtype=object)
print('z-independent (mod p):', np.array_equal(WA % P, WB % P))
print('rank (mod p):', rank_mod(WA))
# grade structure: grade of product basis index (i1,j1,i2,j2) = (i1-j1+3)+(i2-j2+3); word grade a1+k2+a3+k4 ?
grade_pb = {}
for i1 in range(4):
    for j1 in range(4):
        for i2 in range(4):
            for j2 in range(4):
                grade_pb[((i1 * 4 + j1) * 4 + i2) * 4 + j2] = (i1 - j1 + 3) + (i2 - j2 + 3)
ok = True
for n, l in enumerate(labs):
    gs = set(grade_pb[c] for c in range(256) if WA[n, c] % P != 0)
    if gs != {sum(l)}: ok = False
print('each word has single total grade a1+k2+a3+k4:', ok)
# pair-grade structure: does word have definite pair-1 grade?
pg = []
for n, l in enumerate(labs):
    g1 = set((c // 64 - (c // 16) % 4 + 3) for c in range(256) if WA[n, c] % P != 0)
    pg.append((l, sorted(g1)))
print('examples of pair-1 grades per word:', pg[:20])
pickle.dump({'labs': labs, 'W': WA}, open('twopair_W.pkl', 'wb'))
