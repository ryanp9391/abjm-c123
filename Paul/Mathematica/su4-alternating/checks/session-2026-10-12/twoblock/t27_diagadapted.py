import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p
z = z0
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
e1, e4 = np.array(e(1)), np.array(e(4))
def A2c(v): return X.vm(X.vec(v), X.VT)
def F2c(u): return X.vm(X.vec(u), X.Vi)
def show(name, vs):
    c = X.coeffs(X.bra(vs)); print(f'{name:75s}', X.ranks(c))
mz = sum(z) / 4
# unweighted diag-frame twist-adapted: v'_j = p(z_j), u'_j = q(z_j), sum_j p q = 0
show('diag A (1,1,1,1), F q(z_j)=z_j-mean   [block1, block2 = companion vac]', [A2c([1, 1, 1, 1]), F2c([zj - mz for zj in z]), e1, e4])
show('diag A p(z_j)=z_j, F q with sum p q = 0 (q = z^2 - c)        [block1]', [A2c(z), F2c([zj ** 2 - sum(x ** 3 for x in z) / sum(z) for zj in z]), e1, e4])
# Lagrange-weighted = companion flags
for dp in range(4):
    for dq in range(4):
        if dp + dq > 2: continue
        pv = [lam[j] * (1 + sum(Fr(k + 2, 3) * z[j] ** (k + 1) for k in range(dp))) for j in range(4)]
        qv = [(1 + sum(Fr(-k - 1, 5) * z[j] ** (k + 1) for k in range(dq))) for j in range(4)]
        show(f'diag A lam_j p(z_j) deg {dp}, F q(z_j) deg {dq} (auto BPS)  [block1]', [A2c(pv), F2c(qv), e1, e4])
        show(f'   same on both blocks (global)', [A2c(pv), F2c(qv), A2c(pv), F2c(qv)])
