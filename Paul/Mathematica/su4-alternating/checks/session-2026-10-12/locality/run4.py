import sys
from words import *
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4   # sites
a0 = 1/3
ch = build_chain(L)
psi0 = W_state(L, a0)
F = factor_list(ch)
print('nonzero factors', len(F), 'zero dropped', 54 - len(F))
print([f['name'] for f in F if f['s']==0][:30])
# level 1
t = time.time()
B1 = np.stack([apply_f(ch, f, psi0) for f in F], axis=1)
print('level1', B1.shape, time.time() - t)
nz = np.linalg.norm(B1, axis=0)
print('zero bras at level1:', [F[i]['name'] for i in range(len(F)) if nz[i] < 1e-12*np.linalg.norm(psi0)])
S, sv = None, None
def rank_of(M, tol_list=(1e-6,1e-8,1e-10,1e-12)):
    Mn, _ = colnorm(M)
    sv = np.linalg.svd(Mn, compute_uv=False)
    return sv
M = np.concatenate([psi0[:, None], B1], axis=1)
sv = rank_of(M); print('span <=1 factor: sv tail', sv[:3], sv[-5:])
for tol in (1e-6, 1e-8, 1e-10, 1e-12): print(tol, int((sv > tol*sv[0]).sum()))
