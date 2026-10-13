import sys, pickle
from cyc import *
L = int(sys.argv[1]); thp = int(sys.argv[2]); pi = int(sys.argv[3])
cy = Cyclic(L, thp, pi)
Lb = cy.loc_basis_coords()
print('loc dim', Lb.shape[0], flush=True)
ops = np.stack([cy.op(f) for f in cy.F], axis=0)
tg = {k: v for k, v in cy.Tc.items()}
names = [(f['a'], f['s'], f['m']) for f in cy.F]
np.savez_compressed(f'cyc_L{L}_th{thp}_p{pi}.npz', ops=ops, rows=cy.basis.rows, piv=cy.piv, e0=cy.e0, Lb=Lb)
pickle.dump((names, tg, cy.p), open(f'cyc_L{L}_th{thp}_p{pi}.pkl', 'wb'))
print('saved', ops.shape, flush=True)
