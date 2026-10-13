import pickle
from collections import Counter
from fractions import Fraction as Fr
from sov_l4 import WORDS
d = pickle.load(open('gram_exact.pkl', 'rb')); G = d['G']
nz = [(WORDS[i], WORDS[j]) for i in range(256) for j in range(256) if G[i, j] != 0]
c = Counter((sum(x), sum(y)) for x, y in nz)
print('(S_x, S_y) counts:', sorted(c.items()))
# per-site relations: for nonzero entries, pattern of (x_s, y_s) per site
for s in range(4):
    cs = Counter((x[s], y[s]) for x, y in nz)
    print('site', s + 1, sorted(cs.items()))
# per pair (A1,F2) and (A3,F4)
c2 = Counter((x[0] + x[1], y[0] + y[1]) for x, y in nz)
print('pair12 charges', sorted(c2.items()))
c3 = Counter((x[0] + x[3], y[0] + y[3], x[1]+x[2], y[1]+y[2]) for x, y in nz)
print('x(A1)+x(F4) etc', sorted(c3.items())[:20])
# antidiagonal? check G[x, y] with y = (3,3,3,3) - x
anti = sum(1 for x, y in nz if all(x[s] + y[s] == 3 for s in range(4)))
print('# nonzeros with y = 3 - x sitewise:', anti)
diag = sum(1 for x, y in nz if x == y)
print('# nonzeros with y = x:', diag)
