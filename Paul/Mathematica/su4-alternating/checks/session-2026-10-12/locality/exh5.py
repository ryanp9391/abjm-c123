import sys, itertools, time
import numpy as np
from pipeline import *
from span import *
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F0 = factor_list(ch)
F = [f for f in F0 if abs(f['m']) <= 1]
n = len(F)
vecs = [psi0]; labs = [()]
fr = [(psi0, ())]
for l in (1, 2):
    nf = []
    for v, lab in fr:
        last = lab[-1] if lab else 0
        for i in range(last, n):
            nf.append((apply_f(ch, F[i], v), lab + (i,)))
    vecs += [x[0] for x in nf]; labs += [x[1] for x in nf]; fr = nf
W = np.stack(vecs, axis=1); nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, _ = orth_basis(Wn, tol=1e-12)
C = Q.T @ Wn; N = C.shape[1]
print('family', n, 'factors', N, 'words, rank', Q.shape[1], flush=True)
z = np.array([float(x) for x in Z4])
v = np.array([1, 0, 0, a0]); u = np.array([-a0, 0, 0, 1.0]); rest = W_state_rest(L, a0)
T = targets(L)
tg = {"T1 a'=+0.6667": T["T1 a'=+0.6667"], "T2 kappa=+0.3333": T["T2 kappa=+0.3333"],
      'F-site rotation v x uLam': np.kron(np.kron(v, u * z), rest), 'singlet S x w0': np.kron(np.eye(4).ravel(), rest)}
tcs = {k: (Q.T @ t) / np.linalg.norm(Q.T @ t) for k, t in tg.items()}
for kmax in (5,):
    t0 = time.time(); hits = {k: [] for k in tg}
    combos = itertools.combinations(range(N), kmax)
    while True:
        batch = list(itertools.islice(combos, 60000))
        if not batch: break
        idx = np.array(batch)
        M = C[:, idx].transpose(1, 0, 2)          # (b, r, k)
        P = np.linalg.pinv(M, rcond=1e-10)        # (b, k, r)
        for k, tc in tcs.items():
            sol = P @ tc
            res = np.linalg.norm(np.einsum('brk,bk->br', M, sol) - tc, axis=1)
            for j in np.nonzero(res < 1e-7)[0]: hits[k].append((batch[j], res[j]))
    print(f'k = {kmax}: time {time.time()-t0:.0f}s ', {k: len(h) for k, h in hits.items()}, flush=True)
    for k, h in hits.items():
        if h: print('   ', k, [ [' '.join(F[i]['name'] for i in labs[w]) or '1' for w in h[0][0]] ], 'res %.1e' % h[0][1])
