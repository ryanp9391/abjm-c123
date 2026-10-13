import sys, pickle
from pipeline import *
from span import *
L=6
W = np.load('W_L3.npy'); F, labs = pickle.load(open('W_L3.pkl','rb'))
nr = np.linalg.norm(W, axis=0); print('min/max col norms', nr.min(), nr.max())
Wn = W/np.where(nr>0, nr, 1)
t=time.time()
Q, gaps = orth_basis(Wn, tol=1e-13, chunk=1000)
print('rank at tol 1e-13 (chunk-residual svd):', Q.shape[1], time.time()-t)
C = Q.T @ Wn
s = np.linalg.svd(C, compute_uv=False)
print('top', s[:3])
for tol in (1e-4,1e-6,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13):
    print('count s/s0 >', tol, int((s/s[0] > tol).sum()))
i0 = int((s/s[0] > 1e-13).sum())
print('s around:', s[300:i0+3][::10] if i0>300 else s)
np.save('Q_L3_all.npy', Q)
np.save('s_L3.npy', s)
