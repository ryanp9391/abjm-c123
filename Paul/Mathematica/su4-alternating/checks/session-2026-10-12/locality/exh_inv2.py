import sys, itertools, time
import numpy as np
from pipeline import *
from span import *
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
sing = [apply_f(ch, f, psi0) for f in F]
inv = []
for f in F:
    lu, rc = get_inverse_lu(ch, f)
    inv.append(apply_inv(lu, psi0))
vecs = [psi0] + sing + inv
labs = ['1'] + [f['name'] for f in F] + [f['name'] + '^-1' for f in F]
W = np.stack(vecs, axis=1); nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, g = orth_basis(Wn, tol=1e-10)
C = Q.T @ Wn; N = C.shape[1]
print('words', N, 'float rank', Q.shape[1], flush=True)
T = targets(L)
tg = {k: T[k] for k in T if k.startswith('T1') or k.startswith('T2')}
print('residuals outside span:', {k: '%.1e' % resid(Q, v) for k, v in tg.items()})
tcs = {k: (Q.T @ v) / np.linalg.norm(Q.T @ v) for k, v in tg.items()}
# L1 support and OMP on this pool
from sparse import l1, omp
for k in list(tg)[:3] + list(tg)[4:6]:
    sup, x, ok = l1(C, tcs[k], scale_cols=False)
    s2, c2, hist = omp(C, tcs[k], 12)
    print(k, 'L1 support', len(sup) if ok else None, '| OMP residual after 12:', '%.1e' % hist[-1], flush=True)
for kmax in (1, 2, 3, 4):
    t0 = time.time(); hits = {k: [] for k in tg}
    combos = itertools.combinations(range(N), kmax)
    while True:
        batch = list(itertools.islice(combos, 60000))
        if not batch: break
        idx = np.array(batch)
        M = C[:, idx].transpose(1, 0, 2)
        P = np.linalg.pinv(M, rcond=1e-10)
        for k, tc in tcs.items():
            sol = P @ tc
            res = np.linalg.norm(np.einsum('brk,bk->br', M, sol) - tc, axis=1)
            for j in np.nonzero(res < 1e-7)[0]: hits[k].append((batch[j], res[j]))
    print(f'k = {kmax}: {time.time()-t0:.0f}s', {k: len(h) for k, h in hits.items()}, flush=True)
    for k, h in hits.items():
        if h: print('   ', k, [labs[w] for w in h[0][0]], 'res %.1e' % h[0][1])
