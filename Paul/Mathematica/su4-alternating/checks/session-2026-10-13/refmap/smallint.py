"""smallint.py -- exhaustive scan over v, u in {-1,0,1,2}^4 (v.u = 0, both nonzero, up to overall scaling) at (i) z2 = z3 and (ii) generic z: compare with the closed form."""
from helpers import *
import itertools, sys
out = open('smallint_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
vals = [-1, 0, 1, 2]
y = Fr(5, 4)
zz_eq = [Fr(3, 2), y, y, Fr(4, 7)]; zz_gen = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
cands = []
for v in itertools.product(vals, repeat=4):
    if all(x == 0 for x in v): continue
    for u in itertools.product(vals, repeat=4):
        if all(x == 0 for x in u): continue
        if sum(a * b for a, b in zip(v, u)) != 0: continue
        cands.append((tuple(map(Fr, v)), tuple(map(Fr, u))))
P('candidates (v.u=0):', len(cands))
for name, zz in (('z2=z3', zz_eq), ('generic z', zz_gen)):
    ok_list = []; mism = []
    for v, u in cands:
        ok = test(list(v), list(u), zz, kaps=KAPS[:1] + [Fr(1, 3)])[2]
        if name == 'z2=z3':
            if pred(v, u) != ok: mism.append((v, u))
        if ok: ok_list.append((v, u))
    if name == 'z2=z3':
        P(name, ': admissible', len(ok_list), ' mismatches with closed form:', len(mism), mism[:5])
    else:
        # generic z: closed form = blockwise BPS with singleton blocks: v_j u_j = 0 for each j
        blk = [(v, u) for v, u in ok_list if all(a * b == 0 for a, b in zip(v, u))]
        nonblk = [(v, u) for v, u in ok_list if not all(a * b == 0 for a, b in zip(v, u))]
        allblk = [(v, u) for v, u in cands if all(a * b == 0 for a, b in zip(v, u))]
        P(name, ': admissible', len(ok_list), '; of which blockwise (v_j u_j = 0 for all j):', len(blk), 'of', len(allblk), 'blockwise candidates; NON-blockwise admissible:', len(nonblk), [([str(x) for x in v], [str(x) for x in u]) for v, u in nonblk[:10]])
out.close()
