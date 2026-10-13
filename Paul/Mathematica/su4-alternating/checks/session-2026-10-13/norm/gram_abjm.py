"""Python fallback (kernel busy): regularised bra/ket words and Gram matrix at the ABJM point on the eps-chain, exact."""
import sys, time, pickle
SOV = r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/checks/session-2026-10-12/sovformula'
sys.path.insert(0, SOV)
import numpy as np
from fractions import Fraction as Fr
import sov_l4 as S
from regbases import Reg
eps = Fr(sys.argv[1])
phi = [Fr(9, 10), Fr(7, 20), Fr(3, 5), Fr(-11, 10)]; thdir = [Fr(0), Fr(37, 100)]
thF = tuple(eps * t for t in thdir); z = tuple(1 + eps * p for p in phi)
t0 = time.time()
R = Reg(h=Fr(1, 2), thF=thF, z=z)
Wb = R.bra_words(); Yk = R.ket_words()
print('words built %.0f s' % (time.time() - t0), flush=True)
G = S.gram(Wb, Yk)
print('gram %.0f s' % (time.time() - t0), 'nonzeros', sum(1 for x in G.flat if x != 0), flush=True)
ops = dict(bra=[R.pair_ops_bra(0), R.pair_ops_bra(1)], ket=[R.pair_ops_ket(0), R.pair_ops_ket(1)])
pickle.dump(dict(eps=eps, G=G, ops=ops, z=z, thF=thF), open('gram_%s.pkl' % str(eps).replace('/', '_'), 'wb'))
print('saved %.0f s' % (time.time() - t0))
