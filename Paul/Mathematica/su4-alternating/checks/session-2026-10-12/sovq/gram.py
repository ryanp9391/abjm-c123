import pickle, time, itertools
import numpy as np
from fractions import Fraction as Fr
from sov_l4 import *
h, th = H0, TH0
ch = Chain()
t0 = time.time()
KETVAC = (3, 0, 3, 0)
vac = product_vec([unit(i) for i in KETVAC])
def ket_ops(site, lab):
    if site in (0, 2):
        return [ch.t(1, th[site])] * lab
    return [ch.t(lab, th[site] + lab * h)] if lab else []
cols = []
for lab in WORDS:
    v = vac
    for s in (3, 2, 1, 0):
        for M in ket_ops(s, lab[s]):
            v = M.dot(v)
    cols.append(v)
Y = np.array(cols).T
print('ket words built', time.time() - t0)
W = pickle.load(open('words_exact.pkl', 'rb'))['words']
print('rank mod p of ket words', rank_mod_p(Y.T))
t0 = time.time()
# sparse-aware product
Wnz = [[(c, x) for c, x in enumerate(r) if x != 0] for r in W]
G = np.empty((256, 256), dtype=object)
for i in range(256):
    for j in range(256):
        G[i, j] = sum((x * Y[c, j] for c, x in Wnz[i]), Fr(0))
print('Gram built', time.time() - t0)
nz = [(i, j) for i in range(256) for j in range(256) if G[i, j] != 0]
print('nonzeros of G', len(nz))
pickle.dump({'G': G, 'Y': Y}, open('gram_exact.pkl', 'wb'))
