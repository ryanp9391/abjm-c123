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
C = Q.T @ Wn
T = targets(L)
rng = np.random.default_rng(0)
def member(S, tc, tol=1e-9):
    M = C[:, S]
    c, *_ = np.linalg.lstsq(M, tc, rcond=None)
    return np.linalg.norm(M @ c - tc) < tol * np.linalg.norm(tc)
res = {}
for name in ["T1 a'=+0.6667", "T2 kappa=+0.3333"]:
    tc = Q.T @ T[name]; tc /= np.linalg.norm(tc)
    sizes = []; best = None
    for rep in range(60):
        S = list(range(C.shape[1]))
        order = rng.permutation(S)
        for w in order:
            S2 = [s for s in S if s != w]
            if member(S2, tc): S = S2
        sizes.append(len(S))
        if best is None or len(S) < len(best): best = list(S)
    print(name, 'minimal-irredundant support sizes over 60 random orders: min', min(sizes), 'max', max(sizes), 'hist', np.bincount(sizes)[min(sizes):])
    print('  best set:', [label_str(F, labs[s]) for s in best])
    res[name] = best
pickle.dump(res, open('backelim_L2.pkl', 'wb'))
