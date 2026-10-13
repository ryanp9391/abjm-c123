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
r = Q.shape[1]
w0 = W_state_rest(L, a0); R = len(w0); w0n = w0 / (w0 @ w0)
s, U, Vt = localised(Q, w0, None); nloc = int((s > 1 - 1e-8).sum()); Uloc = U[:, :nloc]
def xof(V): return np.einsum('xrn,r->xn', V.reshape(16, R, -1), w0n)
Cn = Q.T @ Wn
P = np.eye(r) - Uloc @ Uloc.T
Nn = P @ Cn; nn = np.linalg.norm(Nn, axis=0)
zero = set(np.nonzero(nn < 1e-8)[0])
active = [j for j in range(W.shape[1]) if j not in zero]
# basis of nonloc coordinates
u2, s2, _ = np.linalg.svd(Nn, full_matrices=False)
Nr = u2[:, :r - nloc].T @ Nn
A = Nr[:, active] / nn[active]
G = A.T @ A
m = len(active)
# individually localised x space
Xloc = np.stack([xof(W[:, [j]])[:, 0] for j in zero], axis=1)
Xloc = Xloc / np.maximum(np.linalg.norm(Xloc, axis=0), 1e-300)
uu, ss, _ = np.linalg.svd(Xloc, full_matrices=False)
Xb = uu[:, :int((ss > 1e-9 * ss[0]).sum())]
print('Xloc rank', Xb.shape[1])
found = []
t0 = time.time()
for i in range(m):
    for j in range(i + 1, m):
        if abs(G[i, j]) > 1 - 1e-8: continue
        q, _ = np.linalg.qr(A[:, [i, j]])
        rem = A[:, j + 1:] - q @ (q.T @ A[:, j + 1:])
        hit = np.nonzero(np.linalg.norm(rem, axis=0) < 1e-8)[0]
        for h in hit:
            l = j + 1 + h
            if abs(G[i, l]) > 1 - 1e-8 or abs(G[j, l]) > 1 - 1e-8: continue
            ids = [active[i], active[j], active[l]]
            M = Nn[:, ids]
            _, sv, Vh = np.linalg.svd(M)
            c = Vh[-1]
            v = W[:, ids] @ c
            x = xof(v[:, None])[:, 0]
            nx = np.linalg.norm(x); nv = np.linalg.norm(v)
            if nx < 1e-9 * np.linalg.norm(W[:, ids] * np.abs(c), axis=0).max(): continue
            xp = x - Xb @ (Xb.T @ x)
            found.append((ids, c, x, np.linalg.norm(xp) / nx))
print('true 3-circuits with nonzero loc part', len(found), time.time() - t0)
good = [f for f in found if f[3] > 1e-8]
print('of which outside Xloc:', len(good))
for f in good[:30]:
    print('  ', ' | '.join(label_str(F, labs[x]) for x in f[0]), ' coeffs', f[1], ' x', {(k//4+1, k%4+1): round(float(f[2][k]), 5) for k in range(16) if abs(f[2][k]) > 1e-10*np.abs(f[2]).max()})
pickle.dump(found, open('circ6_L2.pkl', 'wb'))
