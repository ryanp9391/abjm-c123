import numpy as np
from core import *
R = Ring('float'); ch = Chain(R)
rng = np.random.default_rng(1)
ops = [ch.t(a, u) for a in (1,2,3) for u in (Fr(37,100), Fr(-113,100), Fr(3,2))]
scales = [np.abs(np.linalg.eigvals(o)).max() for o in ops]
best = None
for trial in range(12):
    c = rng.normal(size=len(ops)) + 1j*rng.normal(size=len(ops))
    H = sum(ci/s*o for ci, s, o in zip(c, scales, ops))
    ev = np.linalg.eigvals(H)
    d = np.abs(ev[:,None]-ev[None,:]) + np.eye(256)*1e9
    g = d.min()/np.abs(ev).max()
    print(trial, g)
    if best is None or g > best[0]: best = (g, c)
np.save('Hcoef.npy', best[1]); print('best', best[0])
