import numpy as np
from core import *
from exact import *
R8 = res8(); ps = list(R8.keys())
# B = (v x u) x e1 x e4 : c_x = sum_ij v_i u_j Y[(i,j,0,3), x]
labs = wordsIdx
def idx(i1, i2, i3, i4): return ((i1 * 4 + i2) * 4 + i3) * 4 + i4
rows = [idx(i, j, 0, 3) for i in range(4) for j in range(4)]
Ysub = {p: R8[p][1][rows, :] for p in ps}
nzcols = sorted(set(np.nonzero(Ysub[ps[0]])[1]))
print('labels with nonzero coefficient for (v x u) x vac (generic v,u incl. non-BPS):', len(nzcols))
print([labs[c] for c in nzcols])
# exact values
tab = {}
for c in nzcols:
    for r in range(16):
        vals = {p: int(Ysub[p][r, c]) for p in ps}
        if vals[ps[0]] == 0: continue
        tab[(r, c)] = reconstruct(vals)
fails = [k for k, v in tab.items() if v is None]
print('entries', len(tab), 'reconstruction failures', len(fails))
for c in nzcols:
    ent = [(f'v{r // 4 + 1}u{r % 4 + 1}', tab[(r, c)]) for r in range(16) if (r, c) in tab]
    print(labs[c], ent)
