"""F-only and A-only reduced dual tensors: support and rank-1 (separability) structure per SoV label (exact, mod p, then rationals)."""
import numpy as np, pickle
from core import *
from exact import *
R8 = res8(); ps = list(R8.keys()); p = ps[0]; Y = R8[p][1]
def idx(i1, i2, i3, i4): return ((i1 * 4 + i2) * 4 + i3) * 4 + i4   # 0-based
# F-only: A sites at e1 (index 0), F components i2, i4 in {1,2,3} (e2,e3,e4)  [BPS with e1]
for name, fixed, comps in [('F-only (A=e1; F in span e2,e3,e4)', 'F', [1, 2, 3]), ('A-only (F=e4; A in span e1,e2,e3)', 'A', [0, 1, 2])]:
    labs = set(); T = {}
    for i in comps:
        for j in comps:
            row = idx(0, i, 0, j) if fixed == 'F' else idx(i, 3, j, 3)
            for x in np.nonzero(Y[row])[0]:
                labs.add(x)
    labs = sorted(labs)
    print(name, ': labels', [wordsIdx[x] for x in labs])
    for x in labs:
        Mx = np.array([[Y[idx(0, i, 0, j) if fixed == 'F' else idx(i, 3, j, 3), x] for j in comps] for i in comps], dtype=np.int64)
        print('   label', wordsIdx[x], ' rank of 3x3 (site-s component) matrix:', modrank(Mx, p), ' support', [(comps[a] + 1, comps[b] + 1) for a in range(3) for b in range(3) if Mx[a, b]])
