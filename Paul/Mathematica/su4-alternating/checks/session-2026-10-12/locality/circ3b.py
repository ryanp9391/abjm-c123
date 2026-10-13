import sys, pickle, itertools
from pipeline import *
from span import *
Lsites = int(sys.argv[1])
a0 = 1/3
if Lsites == 6:
    Wfull = np.load('W_L3.npy', mmap_mode='r'); F, labs_all = pickle.load(open('W_L3.pkl', 'rb'))
    W = np.array(Wfull[:, :1225]); labs = labs_all[:1225]
else:
    ch = build_chain(Lsites); psi0 = W_state(Lsites, a0); F = factor_list(ch)
    ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2'), verbose=False)
    W = ws.matrix(); labs = list(ws.labels)
nr = np.linalg.norm(W, axis=0); Wn = W / nr
print('norm range', nr.min(), nr.max())
Q, gaps = orth_basis(Wn, tol=1e-10, chunk=700)
r = Q.shape[1]
print('rank', r)
C = Q.T @ Wn
w0 = W_state_rest(Lsites, a0)
s, U, Vt = localised(Q, w0, None)
print('cos', np.round(s[:8], 8))
nloc = int((s > 1 - 1e-7).sum())
Uloc = U[:, :nloc]
P = np.eye(r) - Uloc @ Uloc.T
Nn = P @ C
nn = np.linalg.norm(Nn, axis=0)
zero = np.nonzero(nn < 1e-8)[0]
print('nloc', nloc, 'individually localised', len(zero), [label_str(F, labs[j]) for j in zero])
u2, s2, _ = np.linalg.svd(Nn, full_matrices=False)
rr = r - nloc
Nr = (u2[:, :rr].T @ Nn)
active = np.nonzero(nn >= 1e-8)[0]
Nu = Nr[:, active] / nn[active]
G = Nu.T @ Nu
iu = np.triu_indices(len(active), 1)
par = np.nonzero(np.abs(G[iu]) > 1 - 1e-8)[0]
print('parallel pairs:', len(par))
for q in par:
    i, j = iu[0][q], iu[1][q]
    print('  ', label_str(F, labs[active[i]]), '|', label_str(F, labs[active[j]]), ' cos %.10f' % G[i, j])
# near-parallel for diagnostics
near = np.nonzero((np.abs(G[iu]) > 1 - 1e-4) & (np.abs(G[iu]) <= 1 - 1e-8))[0]
print('near-parallel (1e-8..1e-4):', len(near))
