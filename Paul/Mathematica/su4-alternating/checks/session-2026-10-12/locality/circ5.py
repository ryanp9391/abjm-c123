import sys, pickle, itertools
from pipeline import *
from span import *
np.set_printoptions(precision=5, suppress=True, linewidth=160)
Lsites = int(sys.argv[1]); a0 = 1/3
if Lsites == 6:
    Wfull = np.load('W_L3.npy', mmap_mode='r'); F, labs_all = pickle.load(open('W_L3.pkl', 'rb'))
    W = np.array(Wfull[:, :1225]); labs = labs_all[:1225]
else:
    ch = build_chain(Lsites); psi0 = W_state(Lsites, a0); F = factor_list(ch)
    ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2'), verbose=False)
    W = ws.matrix(); labs = list(ws.labels)
w0 = W_state_rest(Lsites, a0); R = len(w0); w0n = w0 / (w0 @ w0)
def xof(V): return np.einsum('xrn,r->xn', V.reshape(16, R, -1), w0n)
def idx(s):
    for j, l in enumerate(labs):
        if label_str(F, l) == s: return j
pairs = [('t1(th1-1h)', 't2(th1+0h) t3(th1-1h)'), ('t2(th1-1h)', 't3(th1-1h) t3(th1+0h)'), ('t2(th1+0h)', 't1(th1-1h) t1(th1+0h)'),
         ('t3(th1+0h)', 't1(th1+0h) t2(th1-1h)'), ('t1(th2+1h)', 't2(th2+1h) t3(th2+3h)'), ('t2(th2+1h)', 't1(th2+0h) t1(th2+1h)'),
         ('t2(th2+2h)', 't3(th2+2h) t3(th2+3h)'), ('t3(th2+2h)', 't1(th2+0h) t2(th2+2h)'),
         ('t1(th1-1h) t3(th1+0h)', 't2(th1-1h) t2(th1+0h)'), ('t1(th2+1h) t3(th2+2h)', 't2(th2+1h) t2(th2+2h)')]
for a, b in pairs:
    i, j = idx(a), idx(b)
    # find d minimizing |W_i + d W_j|
    d = -(W[:, i] @ W[:, j]) / (W[:, j] @ W[:, j])
    v = W[:, i] + d * W[:, j]
    nv = np.linalg.norm(v)
    x = xof(v[:, None])[:, 0]
    loc_def = np.linalg.norm(v.reshape(16, R) - np.outer(x, w0))
    print(f'{a:28s} + d*[{b:28s}] d={d: .6e}  |combo|/|W_i|={nv/np.linalg.norm(W[:,i]):.2e}  loc-defect/|combo|={loc_def/max(nv,1e-300):.2e}  |x|={np.linalg.norm(x):.3e}')
