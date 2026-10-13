import sys
from cyc import *
L = int(sys.argv[1]); maxl = int(sys.argv[2])
cy = Cyclic(L, 1, 0)
p = cy.p
Lb = cy.loc_basis_coords()
print('cyclic', cy.r, 'loc', Lb.shape[0], flush=True)
ops = {i: cy.op(f) for i, f in enumerate(cy.F)}
print('ops built', flush=True)
def layered(fidx, maxl):
    B = ModBasis(cy.r, p); B.add(cy.e0[None, :].copy())
    new = cy.e0[None, :].copy()
    res = []
    for l in range(1, maxl + 1):
        cand = np.concatenate([modmat_mul(new, ops[i], p) for i in fidx], axis=0)
        added = []
        for c0 in range(0, cand.shape[0], 400):
            ch = cand[c0:c0 + 400]
            fl = B.add(ch)
            added.append(ch[np.array(fl)])
        new = np.concatenate(added, axis=0)
        mem = {k: B.contains(v) for k, v in cy.Tc.items()}
        res.append((l, B.rank, all(mem[k] for k in mem if k.startswith('T1')), all(mem[k] for k in mem if k.startswith('T2'))))
        if new.shape[0] == 0: break
    return res
fam = {
 'all 48': list(range(len(cy.F))),
 'A-site only (s=1)': [i for i, f in enumerate(cy.F) if f['s'] == 0],
 'F-site only (s=2)': [i for i, f in enumerate(cy.F) if f['s'] == 1],
 't1 only': [i for i, f in enumerate(cy.F) if f['a'] == 1],
 't3 only': [i for i, f in enumerate(cy.F) if f['a'] == 3],
 't2 only': [i for i, f in enumerate(cy.F) if f['a'] == 2],
 't1,t3': [i for i, f in enumerate(cy.F) if f['a'] in (1, 3)],
 't1,t2': [i for i, f in enumerate(cy.F) if f['a'] in (1, 2)],
 't2,t3': [i for i, f in enumerate(cy.F) if f['a'] in (2, 3)],
 '|m|<=2': [i for i, f in enumerate(cy.F) if abs(f['m']) <= 2],
 'm in -1..1': [i for i, f in enumerate(cy.F) if abs(f['m']) <= 1],
 'm in -2..0': [i for i, f in enumerate(cy.F) if -2 <= f['m'] <= 0],
}
for k, fi in fam.items():
    print(k, len(fi), layered(fi, maxl), flush=True)
