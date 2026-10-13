import sys, pickle
from pipeline import *
from span import *
from sparse import *
L = 4; a0 = 1/3
ch = build_chain(L)
psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2','W3AFF','W3AAF'), verbose=False)
Wp = ws.matrix(); labs = list(ws.labels)
nsing = 1 + n
Vinv, Linv, info = build_inverse_words(ch, F, range(n), Wp[:, :nsing], labs[:nsing], verbose=False)
W = np.concatenate([Wp, Vinv], axis=1); labs = labs + Linv
print('words', W.shape)
nr = np.linalg.norm(W, axis=0)
Wn = W / np.where(nr > 0, nr, 1)
Q, gaps = orth_basis(Wn, tol=1e-12)
print('rank', Q.shape[1])
T = targets(L)
w0 = W_state_rest(L, a0)
C = Q.T @ Wn
A = C / 1.0
def show(sup, coef, nrm=nr):
    return ' + '.join('(%.6g)*%s' % (c, label_str(F, labs[j])) for j, c in zip(sup, coef))
for name, tv in T.items():
    tc = Q.T @ tv
    rr = resid(Q, tv)
    print('\n==', name, 'resid outside span %.2e' % rr, flush=True)
    if rr > 1e-8: continue
    # coordinates of unnormalised columns: scale back
    Aun = C * nr[None, :]
    sup, x, ok = l1(Aun, tc)
    if ok:
        print(' L1 support size', len(sup), ' residual', np.linalg.norm(Wn[:, sup] @ (x * nr[sup]) - tv) / np.linalg.norm(tv))
        print('  ', show(sup, x))
    s2, c2, hist = omp(Aun, tc, 12)
    print(' OMP', len(s2), ['%.1e' % h for h in hist])
    print('  ', show(s2, c2))
