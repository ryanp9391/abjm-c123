import sys
from pipeline import *
from span import *
np.set_printoptions(precision=5, suppress=True, linewidth=150)
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2','W3AFF','W3AAF'), verbose=False)
W = ws.matrix(); labs = ws.labels
w0 = W_state_rest(L, a0); R = len(w0); w0n = w0 / (w0 @ w0)
def xof(V): return np.einsum('xrn,r->xn', V.reshape(16, R, -1), w0n)
M = W.reshape(16, R, -1)
inl = xof(W)
d = np.linalg.norm((M - inl[:, None, :] * w0[None, :, None]).reshape(W.shape[0], -1), axis=0) / np.maximum(np.linalg.norm(W, axis=0), 1e-300)
loc = np.nonzero(d < 1e-9)[0]
print(len(loc))
X = inl[:, loc]
for j, k in enumerate(loc):
    print(label_str(F, labs[k]), '  nonzero entries of x[A,F] (A idx,F idx):', {(i//4+1, i%4+1): round(float(X[i, j]), 6) for i in range(16) if abs(X[i, j]) > 1e-10})
# target x
T = targets(L)
for name, tv in list(T.items())[:7]:
    x = xof(tv[:, None])[:, 0]
    sol, res, rk, sv = np.linalg.lstsq(X, x, rcond=None)
    print(name, 'x:', {(i//4+1, i%4+1): round(float(x[i]), 5) for i in range(16) if abs(x[i]) > 1e-12}, 'fit resid', np.linalg.norm(X @ sol - x))
print('rank of X', np.linalg.matrix_rank(X, tol=1e-12))
