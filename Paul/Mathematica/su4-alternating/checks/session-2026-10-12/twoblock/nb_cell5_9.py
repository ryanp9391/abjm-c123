import numpy as np, time
from core import *
from nb_common import *
t0 = time.time()
# ---- cell 5 (exact mod p, two primes) ----
for p in (P1, P2):
    R = Ring('mod', p); ch = Chain(R); th = ch.th; h = ch.h
    def comm(a, b):
        A = ch.t(a, Fr(3, 10)); B = ch.t(b, Fr(-7, 10)); return int(np.count_nonzero(R.red(R.mm(A, B) - R.mm(B, A))))
    print('p', p, 'nonzero entries of commutators (1,3),(1,2),(2,3):', comm(1, 3), comm(1, 2), comm(2, 3))
# float magnitudes for calibration (i) two-block
R = Ring('float'); ch = Chain(R); th = ch.th; h = ch.h
c1 = Fr(1) / (h * (th[0] - th[1]) * (th[0] - th[2] + h) * (th[0] - th[3]))
t1 = ch.t(1, th[0]); Gm = np.linalg.inv(np.array([[float(x) for x in r] for r in ch.G])).T
def calI(vs):
    om = kronvecs(R, vs); lhs = np.linalg.solve(t1.T, om)
    rhs = float(c1) * kronvecs(R, [np.array(vs[0], float) @ Gm] + vs[1:])
    return np.abs(lhs - rhs).max(), np.abs(lhs).max()
vA, uA, vB, uB = [1, 2, -1, 3], [1, 1, 0, -1], [3, 1, 1, 2], [1, -1, 0, -1]
print('c1 =', c1)
print('calib (i) global (float, abs, scale):', calI([vA, uA, vA, uA]))
print('calib (i) two-block (float, abs, scale):', calI([vA, uA, vB, uB]))
Rm = Ring('mod', P1); chm = Chain(Rm)
t1m = chm.t(1, th[0]); t1mi = Rm.inv(t1m); Gmm = Rm.inv(Rm.arr(chm.G)).T.copy()
def calIm(vs):
    om = kronvecs(Rm, vs); lhs = Rm.mm(om[None, :], t1mi)[0]
    rhs = Rm.smul(c1, kronvecs(Rm, [Rm.mm(Rm.arr(vs[0])[None, :], Gmm)[0]] + vs[1:]))
    return int(np.count_nonzero(Rm.red(lhs - rhs)))
print('calib (i) exact mod p: global nonzero residual entries', calIm([vA, uA, vA, uA]), '; two-block', calIm([vA, uA, vB, uB]))
# ---- cell 9 (float64 states, notebook H) ----
z = z0
chd = Chain(R, G=[[z[i] if i == j else 0 for j in range(4)] for i in range(4)])
Vd, gap = states_diag(chd)
V = Rmat(z) @ Vd
M, vac = sov_words(ch)
Wv = wvals(ch, V); Wx = wprod(Wv); ov0 = vac @ V
lhs = M @ V
rel = (np.abs(lhs - Wx * ov0[None, :]).max(axis=0) / np.abs(lhs).max(axis=0)).max()
T13 = ch.t(1, Fr(1, 3)); tv = T13 @ V; tt = tau(T13, V)
eres = (np.abs(tv - V * tt[None, :]).max(axis=0) / np.abs(tv).max(axis=0)).max()
print('cell 9: min gap in sectors %.4g, eig residual t1(1/3) %.2e, max rel residual <x|Psi> vs W<0|Psi> %.2e, min |<0|Psi>|/|Psi| %.4g' % (gap, eres, rel, (np.abs(ov0) / np.linalg.norm(V, axis=0)).min()))
np.savez('nb_states_std.npz', V=V, Vd=Vd, M=M, vac=vac, A1=Wv[0], r2=Wv[1], A3=Wv[2], r4=Wv[3])
print('time', time.time() - t0)
