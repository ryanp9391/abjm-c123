import sys, pickle, itertools
from pipeline import *
from span import *
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2'), verbose=False)
W = ws.matrix(); labs = list(ws.labels)
nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, _ = orth_basis(Wn, tol=1e-12)
r = Q.shape[1]
C = Q.T @ Wn         # r x N coordinates, unit-norm words
w0 = W_state_rest(L, a0)
s, U, Vt = localised(Q, w0, None)
nloc = int((s > 1 - 1e-8).sum())
# orthonormal basis of the localised subspace inside span(Q): columns of U[:, :nloc] (in Q-coordinates)
Uloc = U[:, :nloc]
P = np.eye(r) - Uloc @ Uloc.T
Nn = P @ C            # non-loc part in Q-coordinates (rank r - nloc)
Ll = Uloc.T @ C       # loc coordinates (nloc x N)
print('r', r, 'nloc', nloc)
# individually localised words
nn = np.linalg.norm(Nn, axis=0)
zero = np.nonzero(nn < 1e-9)[0]
print('individually localised:', len(zero))
# basis of nonloc: reduce to rank r-nloc coords
u2, s2, _ = np.linalg.svd(Nn, full_matrices=False)
Nr = (u2[:, :r - nloc].T @ Nn)
Nu = Nr / np.where(nn > 1e-9, nn, 1)   # unit
active = np.nonzero(nn >= 1e-9)[0]
print('active words', len(active))
G = Nu[:, active].T @ Nu[:, active]
iu = np.triu_indices(len(active), 1)
par = np.nonzero(np.abs(G[iu]) > 1 - 1e-9)[0]
print('parallel pairs (k=2 circuits):', len(par))
for q in par[:10]:
    i, j = iu[0][q], iu[1][q]
    print('  ', label_str(F, labs[active[i]]), '|', label_str(F, labs[active[j]]))
# k=3 circuits
A = Nu[:, active]
m = len(active)
circ3 = []
t0 = time.time()
for i in range(m):
    for j in range(i + 1, m):
        if abs(G[i, j]) > 1 - 1e-9: continue
        B = A[:, [i, j]]
        # project all l>j
        q, _ = np.linalg.qr(B)
        rem = A[:, j + 1:] - q @ (q.T @ A[:, j + 1:])
        nrm = np.linalg.norm(rem, axis=0)
        hit = np.nonzero(nrm < 1e-8)[0]
        for h in hit: circ3.append((i, j, j + 1 + h))
print('k=3 circuits', len(circ3), time.time() - t0)
pickle.dump((active, circ3), open('circ3_L2.pkl', 'wb'))
for c in circ3[:20]:
    print('  ', ' | '.join(label_str(F, labs[active[x]]) for x in c))
