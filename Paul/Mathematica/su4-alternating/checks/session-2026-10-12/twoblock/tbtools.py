import numpy as np, pickle
from core import *
_res = pickle.load(open('dual_mod.pkl', 'rb'))

class ModCtx:
    def __init__(self, p=P1):
        self.p = p; self.R = Ring('mod', p)
        self.M, self.Y = _res[p]
        z = z0
        self.V = self.R.arr(Vand(z)); self.Vi = modinv(self.V, p)
        self.VT = self.V.T.copy()
        self.Gc = self.R.arr(Gc_of(z)); self.GcT = self.Gc.T.copy()
        self.Lam = self.R.arr([[z[i] if i == j else 0 for j in range(4)] for i in range(4)])
    def vec(self, v): return self.R.arr(v) if not isinstance(v, np.ndarray) else v % self.p
    def vm(self, v, A): return self.R.mm(self.vec(v)[None, :], A)[0]
    def diag2comp(self, vs):
        """diag-frame site vectors (A,F,A,F) -> companion frame: A: v V^T, F: u V^{-1}"""
        out = []
        for s, v in enumerate(vs):
            out.append(self.vm(v, self.VT if types[s] == 'A' else self.Vi))
        return out
    def bra(self, vs): return kronvecs(self.R, [self.vec(v) for v in vs])
    def coeffs(self, b):
        return self.R.mm(b[None, :] % self.p, self.Y)[0]
    def ranks(self, c):
        p = self.p; c4 = c.reshape(4, 4, 4, 4)
        r12 = modrank(c4.reshape(16, 16), p)
        rs = [modrank(np.moveaxis(c4, s, 0).reshape(4, 64), p) for s in range(4)]
        r13 = modrank(np.transpose(c4, (0, 2, 1, 3)).reshape(16, 16), p)
        return dict(nnz=int(np.count_nonzero(c)), pair=r12, sites=rs, r13_24=r13)
    def r(self, x): return self.R.c(x)
    def rnd(self, rng): return int(rng.integers(1, self.p))
