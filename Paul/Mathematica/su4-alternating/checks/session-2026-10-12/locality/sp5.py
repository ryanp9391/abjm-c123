import sys, pickle, itertools
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
C = Q.T @ Wn
z = np.array([float(x) for x in Z4])
v = np.array([1, 0, 0, a0]); u = np.array([-a0, 0, 0, 1.0])
rest = W_state_rest(L, a0)
tg = {
 'F-site rotation  v x uLam (x) w0': np.kron(np.kron(v, u * z), rest),
 'singlet S (x) w0': np.kron(np.eye(4).ravel(), rest),
 'rho1(Lam^-1)rho2(Lam) psi0': np.kron(np.kron(v / z, u * z), rest),
}
for name, tv in tg.items():
    tc = Q.T @ tv; rr = resid(Q, tv)
    print(name, 'resid outside span', '%.1e' % rr)
    tcn = tc / np.linalg.norm(tc)
    sup, x, ok = l1(C, tcn, scale_cols=False)
    print('   L1 support', len(sup) if ok else None)
    for k in range(1, 9):
        s2, c2, hist = omp(C, tcn, k)
        if hist[-1] < 1e-9:
            print('   OMP exact with', k, 'words:', [label_str(F, labs[j]) for j in s2], np.round(c2, 6)); break
    else:
        print('   OMP residual after 8 words', '%.1e' % hist[-1])
    # exhaustive subsets on singles+identity+inverse? restrict to singles (49)
    for k in range(1, 5):
        hits = subset_search(C[:, :49], tcn, k, maxhits=3, tmax=60)
        if hits and hits[0] != 'timeout':
            print('   exhaustive over {1,singles}: k =', k, 'hit', [label_str(F, labs[j]) for j in hits[0][0]]); break
