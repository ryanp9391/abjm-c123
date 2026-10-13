import numpy as np, time
from core import *
from states import *
t0 = time.time()
R = Ring('float')
z = z0
chd = Chain(R, G=[[z[i] if i == j else 0 for j in range(4)] for i in range(4)])
ch = Chain(R)
evs, Vd = eigenstates_diag(chd)
d = np.abs(evs[:, None] - evs[None, :]) + np.eye(256) * 1e9
print('min gap (rel)', d.min() / np.abs(evs).max())
for a, u in [(1, Fr(1, 3)), (2, Fr(-2, 9)), (3, Fr(5, 13))]:
    tv = chd.t(a, u) @ Vd; tt = tau(chd, a, u, Vd)
    print('diag frame eig residual t_%d' % a, (np.abs(tv - Vd * tt[None, :]).max(axis=0) / np.abs(tv).max(axis=0)).max())
Rm = Rmat(z)
V = Rm @ Vd
for a, u in [(1, Fr(1, 3)), (2, Fr(-2, 9)), (3, Fr(5, 13))]:
    tv = ch.t(a, u) @ V; tt = tau(ch, a, u, V)
    print('companion eig residual t_%d' % a, (np.abs(tv - V * tt[None, :]).max(axis=0) / np.abs(tv).max(axis=0)).max())
M, vac = sov_words(ch)
Ws = wavevals(ch, V); W = Wprod(Ws)
lhs = M @ V; rhs = W * (vac @ V)[None, :]
rel = np.abs(lhs - rhs).max(axis=0) / np.abs(lhs).max(axis=0)
print('max rel residual <x|Psi> = W(x)<0|Psi>, all 256 states:', rel.max())
print('min |<0|Psi>|/|Psi|', (np.abs(vac @ V) / np.linalg.norm(V, axis=0)).min())
print('cond W', np.linalg.cond(W))
np.savez('states_std.npz', V=V, Vd=Vd, evs=evs, WA1=Ws[0], WF2=Ws[1], WA3=Ws[2], WF4=Ws[3], M=M, vac=vac)
print('time', time.time() - t0)
