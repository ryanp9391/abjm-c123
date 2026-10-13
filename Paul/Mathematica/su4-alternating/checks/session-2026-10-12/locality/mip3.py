import sys, time
from cyc import *
import scipy.sparse as sps
from scipy.optimize import milp, LinearConstraint, Bounds
L = int(sys.argv[1]); target = sys.argv[2]; fam = sys.argv[3]; maxl = int(sys.argv[4]); Mbig = float(sys.argv[5]); tl = float(sys.argv[6])
cy = Cyclic(L, 1, 0, verbose=False)
p = cy.p
F = cy.F
if fam == 'all': fidx = list(range(len(F)))
elif fam == 'm1': fidx = [i for i, f in enumerate(F) if abs(f['m']) <= 1]
elif fam == 'm2': fidx = [i for i, f in enumerate(F) if abs(f['m']) <= 2]
ops = {i: cy.op(F[i]) for i in fidx}
# words (<= maxl factors, multisets) as modular coordinate rows
words = [((), cy.e0.copy())]
frontier = [((), cy.e0.copy())]
for l in range(1, maxl + 1):
    newf = []
    for lab, v in frontier:
        last = lab[-1] if lab else -1
        for i in fidx:
            if i < last: continue
            nv = modmat_mul(v[None, :], ops[i], p)[0]
            newf.append((lab + (i,), nv))
    words += newf; frontier = newf
print('words', len(words), flush=True)
# dedupe words that are the zero vector
Wm = np.stack([w[1] for w in words], axis=1)   # r x N modular
# use REAL-valued representation for MIP: convert modular to float is meaningless -> need exact rational/real coordinates.
