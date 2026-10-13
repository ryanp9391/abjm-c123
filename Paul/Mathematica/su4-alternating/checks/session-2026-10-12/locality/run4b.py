import sys
from pipeline import *
from span import *
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
classes = tuple(sys.argv[2].split(',')) if len(sys.argv) > 2 else ('W1', 'W2', 'W3AFF', 'W3AAF')
a0 = 1/3
ch = build_chain(L)
psi0 = W_state(L, a0)
F = factor_list(ch)
t = time.time()
ws, B1, B2 = build_positive_words(ch, psi0, F, classes)
W = ws.matrix()
print('N words', W.shape, 'build time', time.time() - t, flush=True)
Q, gaps = orth_basis(W)
print('rank', Q.shape[1])
allg = np.sort(np.concatenate(gaps))[::-1]
print('sv of residuals around the cut:', allg[allg > 1e-14][-5:], 'count >1e-9', (allg > 1e-9).sum())
w0 = W_state_rest(L, a0)
s, U, Vt = localised(Q, w0, None)
print('cos angles to localised space:', np.round(s, 6))
for k, tv in targets(L).items():
    print(k, 'resid', resid(Q, tv))
