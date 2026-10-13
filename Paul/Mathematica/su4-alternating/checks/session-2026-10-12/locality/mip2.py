import sys, pickle, itertools, time
from pipeline import *
from span import *
from scipy.optimize import milp, LinearConstraint, Bounds
L = 4; a0 = 1/3
ch = build_chain(L); psi0 = W_state(L, a0)
F = factor_list(ch); n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2'), verbose=False)
W = ws.matrix(); labs = list(ws.labels)
nr = np.linalg.norm(W, axis=0); Wn = W / nr
Q, _ = orth_basis(Wn, tol=1e-12)
C = Q.T @ Wn
N = C.shape[1]; r = C.shape[0]
T = targets(L)
name = sys.argv[1]; Mbig = float(sys.argv[2]); tl = float(sys.argv[3])
tc = Q.T @ T[name]; tc /= np.linalg.norm(tc)
# variables: c (N), z (N); minimize sum z ; C c = tc ; -M z <= c <= M z
cost = np.concatenate([np.zeros(N), np.ones(N)])
A_eq = np.concatenate([C, np.zeros((r, N))], axis=1)
I = np.eye(N)
A_ub = np.concatenate([np.concatenate([I, -Mbig * I], axis=1), np.concatenate([-I, -Mbig * I], axis=1)], axis=0)
cons = [LinearConstraint(A_eq, tc, tc), LinearConstraint(A_ub, -np.inf, 0)]
integrality = np.concatenate([np.zeros(N), np.ones(N)])
bounds = Bounds(np.concatenate([-Mbig * np.ones(N), np.zeros(N)]), np.concatenate([Mbig * np.ones(N), np.ones(N)]))
t0 = time.time()
res = milp(cost, constraints=cons, integrality=integrality, bounds=bounds, options=dict(time_limit=tl, disp=False, mip_rel_gap=0.0))
print(name, 'status', res.status, res.message, 'time', time.time() - t0)
if res.x is not None:
    c = res.x[:N]; z = res.x[N:]
    sup = np.nonzero(np.abs(c) > 1e-9)[0]
    print('support', len(sup), 'obj', res.fun, 'resid', np.linalg.norm(C @ c - tc))
    for j in sup: print('   %12.5g  %s' % (c[j] / nr[j], label_str(F, labs[j])))
