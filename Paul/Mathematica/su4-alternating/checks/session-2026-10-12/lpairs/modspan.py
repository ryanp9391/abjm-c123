"""incremental row-echelon basis mod P for Krylov closures"""
import numpy as np
from eng import P

class Span:
    def __init__(self, D):
        self.rows = []      # reduced rows (pivot normalised to 1)
        self.piv = []
        self.D = D
    def reduce(self, v):
        v = v.reshape(-1) % P
        v = v.copy()
        for r, c in zip(self.rows, self.piv):
            if v[c]:
                v = (v - (v[c] * r) % P) % P
        return v
    def add(self, v):
        w = self.reduce(v)
        nz = np.nonzero(w)[0]
        if len(nz) == 0: return False
        c = nz[0]; w = (w * pow(int(w[c]), P - 2, P)) % P
        # keep fully reduced: eliminate c from existing rows
        for i, r in enumerate(self.rows):
            if r[c]: self.rows[i] = (r - (r[c] * w) % P) % P
        self.rows.append(w); self.piv.append(c)
        return True
    def contains(self, v):
        return not self.reduce(v).any()
    def dim(self): return len(self.rows)
