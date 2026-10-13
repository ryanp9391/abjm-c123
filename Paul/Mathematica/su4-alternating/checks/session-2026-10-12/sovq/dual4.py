import pickle, itertools
import mpmath as mp
import numpy as np
from fractions import Fraction as Fr
from qfun import *
from poly import *
from sov_l4 import WORDS, vandermonde, minv, mT, Z0
mp.mp.dps = 50
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
ms = pickle.load(open('measure_exact.pkl', 'rb')); M = ms['M']; Gb = ms['Gb']
we = pickle.load(open('words_exact.pkl', 'rb')); W = we['words']
ge = pickle.load(open('gram_exact.pkl', 'rb')); Y = ge['Y']
tomp = lambda x: mp.mpf(x.numerator) / x.denominator
# diagonal-frame eigenvectors -> full 256 vectors; R maps diag kets to companion kets
V = vandermonde(Z0); VmT = mT(minv(V))
def kron(ms_):
    out = np.array([[Fr(1)]], dtype=object)
    for m in ms_: out = np.kron(out, np.array(m, dtype=object))
    return out
R = kron([VmT, V, VmT, V])
Rmp = mp.matrix([[tomp(R[i, j]) for j in range(256)] for i in range(256)])
# left eigenvectors in the diag frame: per block, inverse of the right eigenvector matrix
blocks = data['blocks']
right = [None] * 256; left = [None] * 256
byblock = {}
for n, d in enumerate(data['states']): byblock.setdefault(tuple(d['idx']), []).append(n)
for idx, ns in byblock.items():
    Rb = mp.matrix([[data['states'][n]['vec'][r] for n in ns] for r in range(len(idx))])
    Lb = mp.inverse(Rb)
    for k, n in enumerate(ns):
        v = [mp.mpf(0)] * 256; l = [mp.mpf(0)] * 256
        for r, c in enumerate(idx):
            v[c] = Rb[r, k]; l[c] = Lb[k, r]
        right[n] = v; left[n] = l
Rinv = mp.inverse(Rmp)
Wm = mp.matrix([[tomp(x) for x in r] for r in W])
Ym = mp.matrix([[tomp(Y[i, j]) for j in range(256)] for i in range(256)])
vac_bra = Wm[0, :]   # word (0,0,0,0) = <0|
ketvac = Ym[:, 0]    # ket word (0,0,0,0) = |0bar>
sample = list(range(0, 256, 17))
mx1 = mx2 = 0
Bvals = {}
for n in sample:
    psi = Rmp * mp.matrix(right[n])            # companion ket
    psiL = mp.matrix(left[n]).T * Rinv          # companion bra (row)
    ov = Wm * psi                                # <x|Psi>
    ovL = psiL * Ym                              # <PsiL|y>
    p0 = ov[0]; pL0 = ovL[0]
    # completeness: <PsiL|Psi> = sum_{x,y} <PsiL|y> M_{y,x} <x|Psi>   (M_{y,x} with y = ket label)
    tot = mp.mpf(0)
    for (yb, x), val in M.items():   # stored as M[(xp, x)] with ket label y = 3 - xp
        y = tuple(3 - v for v in yb)
        tot += ovL[WORDS.index(y)] * tomp(val) * ov[WORDS.index(x)]
    norm = (psiL * psi)[0]
    mx1 = max(mx1, abs(tot / norm - 1))
    Bvals[n] = norm / (p0 * pL0)
print('completeness with exact measure: max |sum/<PsiL|Psi> - 1| =', mp.nstr(mx1, 3))
pickle.dump({'right': right, 'left': left}, open('lr_vecs.pkl', 'wb'))
