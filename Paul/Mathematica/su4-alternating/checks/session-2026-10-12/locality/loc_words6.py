import pickle, numpy as np
from span import W_state_rest
from words import *
from pipeline import label_str
W = np.load('W_L3.npy'); F, labs = pickle.load(open('W_L3.pkl', 'rb'))
a0 = 1/3
w0 = W_state_rest(6, a0); R = len(w0); wh = w0 / np.linalg.norm(w0)
res = []
for c0 in range(0, W.shape[1], 2000):
    V = W[:, c0:c0 + 2000]
    M = V.reshape(16, R, -1)
    inl = np.einsum('xrn,r->xn', M, wh)
    nt = np.linalg.norm(V, axis=0)
    out = np.linalg.norm((M - inl[:, None, :] * wh[None, :, None]).reshape(V.shape[0], -1), axis=0)
    for j in range(V.shape[1]):
        if nt[j] > 0 and out[j] / nt[j] < 1e-9:
            res.append(c0 + j)
print('individually localised (<=3 factors, task classes) at L=3:', len(res))
for j in res: print('  ', label_str(F, labs[j]))
