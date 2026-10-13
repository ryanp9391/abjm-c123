# relation between the 4-site SoV words and the tensor product of the two local 2-site SoV bases
import numpy as np, pickle
from fractions import Fraction
from itertools import combinations
from fastchain import *
d = pickle.load(open('twopair_W.pkl', 'rb')); labs = d['labs']; W4 = d['W']
hN = Fraction(1, 2); thN = [Fraction(2, 7), Fraction(-3, 5), Fraction(8, 11), Fraction(-11, 8)]
zA = [Fraction(3, 2), Fraction(5, 4), Fraction(2, 3), Fraction(4, 7)]
cv = [sum(np.prod(s) for s in combinations(zA, k)) for k in range(1, 5)]
Rm = Ring('mod')
def W2(thpair):
    ch = FastChain(['A', 'F'], thpair, hN, cv, Rm)
    vac = Rm.zeros(16); vac[3] = 1
    rows = []
    for a in range(4):
        b = ch.bra_t(vac, a, Rm.red(ch.th[0] - ch.h))
        for k in range(4):
            rows.append(b); b = ch.bra_t(b, 1, ch.th[1])
    return np.array(rows, dtype=object)
W12 = W2(thN[:2]); W34 = W2(thN[2:])
Wloc = np.kron(W12, W34) % P     # row index (a1 k2)(a3 k4) matches labs order (sorted tuples) -> check
lab_loc = [(a1, k2, a3, k4) for a1 in range(4) for k2 in range(4) for a3 in range(4) for k4 in range(4)]
assert lab_loc == labs
Tm = (W4.dot(inv_mod(Wloc))) % P      # <x|_4 = sum_y Tm[x,y] <y|_loc
nz = Tm % P != 0
print('nonzeros of transfer matrix T (256x256):', int(nz.sum()))
# structure: is T upper/lower triangular w.r.t. some natural order?
idx = {l: n for n, l in enumerate(labs)}
def tri_check(key):
    # T[x,y] != 0 only if key(y) >= key(x)?  report counts of violations both ways
    up = sum(1 for x in range(256) for y in range(256) if nz[x, y] and key(labs[y]) < key(labs[x]))
    lo = sum(1 for x in range(256) for y in range(256) if nz[x, y] and key(labs[y]) > key(labs[x]))
    return up, lo
for name, key in [('pair-2 label a3+k4', lambda l: l[2] + l[3]), ('pair-1 label a1+k2', lambda l: l[0] + l[1]),
                  ('total grade', lambda l: sum(l)), ('a3', lambda l: l[2]), ('k4', lambda l: l[3]), ('a1', lambda l: l[0]), ('k2', lambda l: l[1])]:
    print(name, ': (#entries with key(y)<key(x), #with key(y)>key(x)) =', tri_check(key))
# diagonal
print('diagonal entries nonzero:', all(nz[x, x] for x in range(256)))
# within total grade: T maps word x to local words y with same total grade?
print('total grade conserved by T:', tri_check(lambda l: sum(l)) == (0, 0))
pickle.dump({'T': Tm, 'W12': W12, 'W34': W34}, open('twopair_T.pkl', 'wb'))
# pair-1 grade transfer: for each x, the set of (pair-1 grade of y) - (pair-1 grade of x)
from collections import Counter
cnt = Counter()
for x in range(256):
    for y in range(256):
        if nz[x, y]: cnt[(labs[y][0] + labs[y][1]) - (labs[x][0] + labs[x][1])] += 1
print('distribution of (pair-1 grade shift y - x):', sorted(cnt.items()))
