import numpy as np, time, pickle
from core import *
t0 = time.time()
res = {}
for p in (P1, P2):
    R = Ring('mod', p); ch = Chain(R)
    M, vac = sov_words(ch)
    Y = R.inv(M)
    assert np.all(R.mm(M, Y) == np.eye(256, dtype=np.int64))
    res[p] = (M, Y)
    print(p, 'inverse ok', time.time() - t0)
pickle.dump(res, open('dual_mod.pkl', 'wb'))
M, Y = res[P1]; M2, Y2 = res[P2]
nzY = (Y != 0); nzY2 = (Y2 != 0)
print('Y nonzero pattern same for both primes:', np.array_equal(nzY, nzY2), ' # nonzero', nzY.sum())
print('M nonzero', (M != 0).sum(), (M2 != 0).sum())
# pattern: rows of Y = product-basis index (i1 i2 i3 i4), columns = labels (a1 k2 a3 k4)
lab = wordsIdx
cnt_per_label = nzY.sum(axis=0)
print('nonzero components per label: min/max', cnt_per_label.min(), cnt_per_label.max())
# block structure: for each label, which (i1,i2) and (i3,i4) appear
Y4 = nzY.reshape(4, 4, 4, 4, 4, 4, 4, 4)  # i1 i2 i3 i4 a1 k2 a3 k4
# marginal: does (i1,i2) support depend only on (a1,k2)?
S12 = Y4.any(axis=(2, 3))  # i1 i2 a1 k2 a3 k4
S34 = Y4.any(axis=(0, 1))  # i3 i4 a1 k2 a3 k4
dep12 = all(np.array_equal(S12[:, :, a1, k2, a3, k4], S12[:, :, a1, k2, 0, 0]) for a1 in range(4) for k2 in range(4) for a3 in range(4) for k4 in range(4))
print('(i1,i2)-support depends only on (a1,k2):', dep12)
for a1 in range(4):
    for k2 in range(4):
        sup = [(i + 1, j + 1) for i in range(4) for j in range(4) if S12[i, j, a1, k2].any()]
        print('label (a1,k2)=', (a1, k2), ' pair-1 components (union over a3,k4):', len(sup))
