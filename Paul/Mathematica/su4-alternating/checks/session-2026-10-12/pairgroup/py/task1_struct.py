# structure of <0| t_k(theta_s + m h) on the companion vacuum: which site-sets get excited
import numpy as np, itertools, sys
from collections import defaultdict
from pgbase import *
h = Fr(1,2); th = [Fr(2,7), Fr(-3,5), Fr(8,11), Fr(-11,8)]
if len(sys.argv) > 1 and sys.argv[1] == 'abjm':
    th = [Fr(-3,5)+1, Fr(-3,5), Fr(8,11)+1, Fr(8,11)]
z = [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
Gc = companion(z)
ch = Chain(['A','F','A','F'], th, h, Gc)
vidx = (0,3,0,3)
ref = np.zeros((4,)*4, dtype=object); ref[:] = Fr(0); ref[vidx] = Fr(1)
for k in (1,2,3):
    for s in range(2):
        for m in range(-3, 5):
            u0 = th[s] + m*h
            out = ch.apply_t(k, u0, ref, reg=True)
            sets = defaultdict(int)
            for idx in itertools.product(range(4), repeat=4):
                if out[idx] != 0:
                    S = tuple(i+1 for i in range(4) if idx[i] != vidx[i])
                    sets[S] += 1
            if not sets: print('t%d(th%d%+dh): ZERO' % (k, s+1, m)); continue
            print('t%d(th%d%+dh):' % (k, s+1, m), dict(sets))
