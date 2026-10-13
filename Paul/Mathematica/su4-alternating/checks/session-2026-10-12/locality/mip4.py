import sys, time, itertools
from pipeline import *
from span import *
import scipy.sparse as sps
from scipy.optimize import milp, LinearConstraint, Bounds
target = sys.argv[1]; mmax = int(sys.argv[2]); Mbig = float(sys.argv[3]); tl = float(sys.argv[4]); maxl = int(sys.argv[5]) if len(sys.argv) > 5 else 2
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F0 = factor_list(ch)
F = [f for f in F0 if abs(f['m']) <= mmax]
n = len(F)
# words up to maxl factors
vecs = [psi0]; labs = [()]
frontier = [(psi0, ())]
for l in range(1, maxl + 1):
    nf = []
    for v, lab in frontier:
        last = lab[-1] if lab else -1
        for i in range(last if lab else 0, n):
            if lab and i < last: continue
            nf.append((apply_f(ch, F[i], v), lab + (i,)))
    vecs += [x[0] for x in nf]; labs += [x[1] for x in nf]; frontier = nf
W = np.stack(vecs, axis=1)
nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, _ = orth_basis(Wn, tol=1e-12)
C = Q.T @ Wn
N = C.shape[1]; r = C.shape[0]
print('family factors', n, 'words', N, 'rank', r, flush=True)
T = targets(L)
tc = Q.T @ T[target]; res0 = resid(Q, T[target]); print('resid outside span', res0)
tc /= np.linalg.norm(tc)
cost = np.concatenate([np.zeros(N), np.ones(N)])
A_eq = sps.csr_matrix(np.concatenate([C, np.zeros((r, N))], axis=1))
I = sps.identity(N, format='csr')
A_ub = sps.vstack([sps.hstack([I, -Mbig * I]), sps.hstack([-I, -Mbig * I])]).tocsr()
cons = [LinearConstraint(A_eq, tc, tc), LinearConstraint(A_ub, -np.inf, 0)]
integrality = np.concatenate([np.zeros(N), np.ones(N)])
bounds = Bounds(np.concatenate([-Mbig * np.ones(N), np.zeros(N)]), np.concatenate([Mbig * np.ones(N), np.ones(N)]))
t0 = time.time()
res = milp(cost, constraints=cons, integrality=integrality, bounds=bounds, options=dict(time_limit=tl, disp=False, mip_rel_gap=0.0))
print(target, 'status', res.status, res.message, 'time %.0f' % (time.time() - t0))
if res.x is not None:
    c = res.x[:N]
    sup = np.nonzero(np.abs(c) > 1e-9)[0]
    print('support', len(sup), 'obj', res.fun, 'resid', np.linalg.norm(C @ c - tc))
    for j in sup: print('   %12.5g  %s' % (c[j] / nr[j], ' '.join(F[i]['name'] for i in labs[j]) or '1'))
