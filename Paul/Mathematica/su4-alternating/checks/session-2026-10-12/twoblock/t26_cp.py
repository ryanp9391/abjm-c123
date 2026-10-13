"""CP rank (minimal # of products of single-site SoV integrals) of the physical coefficient tensors c_{a1 k2 k4} (a3=0)."""
import numpy as np
from core import *
d = np.load('states_std.npz'); M = d['M']
Yf = np.linalg.inv(M)
z = z0
lam = np.array([float(Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j])) for j in range(4)])
Vf = np.array([[float(x) for x in row] for row in Vand(z)]); Vi = np.linalg.inv(Vf)
vz4 = Vf[:, 3]; l1 = Vi[0]; e1, e4 = np.eye(4)[0], np.eye(4)[3]
pair2 = np.kron(e1, e4)
pieces = {
  'b^1 piece': np.kron(np.kron(lam[0] * vz4, e4) + np.kron(e1, -l1), pair2),
  'b^2 piece': np.kron(np.kron(lam[0] * vz4, -l1), pair2),
}
def als(T, r, iters=3000, seeds=20):
    best = np.inf
    n = T.shape
    for s in range(seeds):
        rng = np.random.default_rng(s)
        A = [rng.normal(size=(n[i], r)) + 1j * rng.normal(size=(n[i], r)) for i in range(3)]
        for it in range(iters):
            for m in range(3):
                o = [i for i in range(3) if i != m]
                KR = np.einsum('ir,jr->ijr', A[o[0]], A[o[1]]).reshape(-1, r)
                Tm = np.moveaxis(T, m, 0).reshape(n[m], -1)
                A[m] = np.linalg.lstsq(KR, Tm.T, rcond=None)[0].T
        rec = np.einsum('ir,jr,kr->ijk', *A)
        res = np.linalg.norm(rec - T) / np.linalg.norm(T)
        best = min(best, res)
        if best < 1e-10: break
    return best
for name, B in pieces.items():
    c = B @ Yf
    c4 = c.reshape(4, 4, 4, 4)
    print(name, 'max |c| outside a3=0, k4<=2:', np.abs(c4[:, :, 1:, :]).max() / np.abs(c).max(), np.abs(c4[:, :, :, 3]).max() / np.abs(c).max())
    T = c4[:, :, 0, :3].copy()
    # normalise per-site scales (diagonal rescalings do not change CP rank)
    for ax in range(3):
        sc = np.sqrt((np.abs(np.moveaxis(T, ax, 0)) ** 2).reshape(T.shape[ax], -1).sum(axis=1))
        sh = [1, 1, 1]; sh[ax] = -1
        T = T / sc.reshape(sh)
    for r in range(3, 8):
        print('   CP rank test r=%d: best relative residual %.2e' % (r, als(T, r)))
