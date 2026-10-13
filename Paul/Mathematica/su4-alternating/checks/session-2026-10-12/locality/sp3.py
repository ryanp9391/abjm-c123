import sys, pickle
from pipeline import *
from span import *
from sparse import *
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2'), verbose=False)
W = ws.matrix(); labs = list(ws.labels)
nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, _ = orth_basis(Wn, tol=1e-12)
r = Q.shape[1]; print('rank', r, 'words', W.shape[1])
C = Q.T @ Wn
T = targets(L)
def irl1(C, tc, iters=15, eps=1e-3):
    N = C.shape[1]
    wgt = np.ones(N)
    best = None
    for it in range(iters):
        # min sum w_j |x_j|  s.t. C x = tc  -> substitute y = w x
        Cs = C / wgt
        sup, x, ok = l1(Cs, tc, scale_cols=False)
        if not ok: return None
        xx = np.zeros(N); xx[sup] = x / wgt[sup]
        wgt = 1.0 / (np.abs(xx) + eps)
        if best is None or len(sup) < len(best[0]): best = (sup, xx[sup])
        eps = max(eps * 0.5, 1e-9)
    return best
for name in ["T1 a'=-0.6667", "T1 a'=+0.6667", "T2 kappa=+0.3333", "T2 kappa=-0.4"]:
    tv = T[name]; tc = Q.T @ tv
    tc = tc / np.linalg.norm(tc)
    sup, x, ok = l1(C, tc, scale_cols=False)
    print(name, 'L1 support', len(sup), 'resid', np.linalg.norm(C[:, sup] @ x - tc))
    b = irl1(C, tc)
    print('   reweighted L1 support', len(b[0]), 'resid', np.linalg.norm(C[:, b[0]] @ b[1] - tc))
    for k in range(1, 8):
        s2, c2, hist = omp(C, tc, k)
        if hist[-1] < 1e-9: print('   OMP exact with', k); break
    else:
        print('   OMP no exact up to 7, hist', ['%.1e' % h for h in hist])
