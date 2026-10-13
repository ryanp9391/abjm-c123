import sys
from pipeline import *
from span import *
L = int(sys.argv[1]); 
a0 = 1/3
ch = build_chain(L)
psi0 = W_state(L, a0)
F = factor_list(ch)
n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1','W2','W3AFF','W3AAF'), verbose=False)
W = ws.matrix(); labs = ws.labels
w0 = W_state_rest(L, a0); wh = w0/np.linalg.norm(w0)
R = len(w0)
def loc_defect(V):
    # V (dim, N): relative norm of component outside C^16 x w0
    M = V.reshape(16, R, -1)
    inl = np.einsum('xrn,r->xn', M, wh)
    nt = np.linalg.norm(V, axis=0)
    out = np.linalg.norm((M - inl[:, None, :] * wh[None, :, None]).reshape(V.shape[0], -1), axis=0)
    return out / np.where(nt>0, nt, 1), nt
d, nt = loc_defect(W)
loc = np.nonzero((d < 1e-9) & (nt > 1e-12*np.linalg.norm(psi0)))[0]
print('N words', W.shape[1], 'individually localised (nonzero):', len(loc))
for j in loc[:60]:
    print(label_str(F, labs[j]), ' |bra|=%.3e' % nt[j])
by = {}
for j in loc: by.setdefault(len(labs[j]), []).append(j)
print({k: len(v) for k, v in by.items()})
