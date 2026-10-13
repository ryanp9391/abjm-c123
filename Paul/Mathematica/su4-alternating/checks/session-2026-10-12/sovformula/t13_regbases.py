import time, pickle
from regbases import *
T0 = time.time()
R = Reg()
Wb = R.bra_words(); print('bra words built %.0f s' % (time.time() - T0), flush=True)
print('rank regularised bra basis (mod p):', S.rank_mod_p(Wb), flush=True)
Yk = R.ket_words(); print('ket words built %.0f s' % (time.time() - T0), flush=True)
print('rank regularised ket basis (mod p):', S.rank_mod_p(Yk.T), flush=True)
G = S.gram(Wb, Yk)
nz = [(S.WORDS[i], S.WORDS[j]) for i in range(256) for j in range(256) if G[i, j] != 0]
print('Gram nonzeros:', len(nz))
from collections import Counter
print('charge sums S_x + S_y of nonzeros:', Counter(sum(x) + sum(y) for x, y in nz))
# pair factorisation test: is G = G_pair1 (x) G_pair2 in the (a,k)x(m,n) labels?
idx = {w: i for i, w in enumerate(S.WORDS)}
pickle.dump(dict(Wb=Wb, Yk=Yk, G=G), open('regbases_L2.pkl', 'wb'))
print('time %.0f s' % (time.time() - T0))
