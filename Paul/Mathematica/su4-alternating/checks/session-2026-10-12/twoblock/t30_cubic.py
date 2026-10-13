import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; R = X.R
z = z0
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
I = np.eye(4, dtype=np.int64)
def polyG(A, roots):
    out = I.copy()
    for r in roots: out = R.mm(out, (A - X.r(r) * I) % p)
    return out
vz4 = X.vm(e(4), X.VT); l1 = X.vm(e(1), X.Vi)
print('v_z4 == e1 prod_{m!=4}(G^T - z_m):', np.array_equal(vz4, X.vm(e(1), polyG(X.GcT, [z[0], z[1], z[2]]))))
print('l1 == lam1 e4 prod_{m!=1}(G - z_m):', np.array_equal(l1, (X.r(lam[0]) * X.vm(e(4), polyG(X.Gc, [z[1], z[2], z[3]]))) % p))
