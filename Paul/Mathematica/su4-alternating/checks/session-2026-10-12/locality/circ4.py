import sys, pickle, itertools
from pipeline import *
from span import *
np.set_printoptions(precision=5, suppress=True, linewidth=160)
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2'), verbose=False)
W = ws.matrix(); labs = list(ws.labels)
nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, _ = orth_basis(Wn, tol=1e-12)
r = Q.shape[1]; C = Q.T @ Wn
w0 = W_state_rest(L, a0); R = len(w0); w0n = w0 / (w0 @ w0)
s, U, Vt = localised(Q, w0, None); nloc = int((s > 1 - 1e-8).sum()); Uloc = U[:, :nloc]
P = np.eye(r) - Uloc @ Uloc.T
Nn = P @ C; nn = np.linalg.norm(Nn, axis=0)
def xof(V): return np.einsum('xrn,r->xn', V.reshape(16, R, -1), w0n)
gens = []   # (cost, description, x vector(16), words with coefficients)
zero = np.nonzero(nn < 1e-8)[0]
for j in zero:
    x = xof(W[:, [j]])[:, 0]
    if np.linalg.norm(x) > 1e-14: gens.append((1, [(labs[j], 1.0)], x))
# parallel pairs
active = np.nonzero(nn >= 1e-8)[0]
u2, s2, _ = np.linalg.svd(Nn, full_matrices=False)
Nr = u2[:, :r - nloc].T @ Nn
Nu = Nr[:, active] / nn[active]
G = Nu.T @ Nu
iu = np.triu_indices(len(active), 1)
par = np.nonzero(np.abs(G[iu]) > 1 - 1e-8)[0]
for q in par:
    i, j = iu[0][q], iu[1][q]
    wi, wj = active[i], active[j]
    # combo W_i*c + W_j*d with nonloc cancel: c*Nn_i + d*Nn_j = 0 -> c = 1, d = -(Nn_i.Nn_j)/(Nn_j.Nn_j)
    d = -(Nn[:, wi] @ Nn[:, wj]) / (Nn[:, wj] @ Nn[:, wj])
    v = W[:, wi] + d * W[:, wj]
    x = xof(v[:, None])[:, 0]
    # check localised
    out = np.linalg.norm(v.reshape(16, R) - np.outer(x, w0))
    gens.append((2, [(labs[wi], 1.0), (labs[wj], d)], x))
    print('circuit', label_str(F, labs[wi]), '|', label_str(F, labs[wj]), 'd=%.5f' % d, 'loc-defect %.1e' % (out / np.linalg.norm(v)), 'x nonzero:', {(k//4+1, k%4+1): round(float(x[k]),5) for k in range(16) if abs(x[k])>1e-12*np.abs(x).max()})
print('generators', len(gens))
X = np.stack([g[2] for g in gens], axis=1)
print('rank of generator x-space', np.linalg.matrix_rank(X, tol=1e-12*np.abs(X).max()))
T = targets(L)
cost = np.array([g[0] for g in gens])
for name, tv in T.items():
    xt = xof(tv[:, None])[:, 0]
    best = None
    for k in range(1, 7):
        for sub in itertools.combinations(range(len(gens)), k):
            M = X[:, sub]
            c, *_ = np.linalg.lstsq(M, xt, rcond=None)
            if np.linalg.norm(M @ c - xt) < 1e-9 * np.linalg.norm(xt):
                tot = cost[list(sub)].sum()
                if best is None or tot < best[0]: best = (tot, sub, c)
        if best is not None and k >= 1 and best[0] <= k: break
    if best:
        print(name, 'min cost (words):', best[0], 'generators:', [(' + '.join('%s*(%.4g)' % (label_str(F, l), c) for l, c in gens[i][2-1] ) if False else [label_str(F, l) for l, c in gens[i][1]]) for i in best[1]])
    else:
        print(name, 'not in span of generators with <= 6 generators')
