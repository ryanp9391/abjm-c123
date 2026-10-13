import time, itertools, pickle
import numpy as np
from fractions import Fraction as Fr
from sov_l4 import *

h, th, z = H0, TH0, Z0
t0 = time.time()
ch = Chain()  # exact companion
T1 = ch.t(1, th[0])
print('exact t1 build', time.time() - t0)

# (i) <Omega| t1(th1)^{-1} = c1 <Omega| rho1(G^{-T})  <=>  <Omega| = c1 <Omega| rho1(G^{-T}) t1(th1)
G = Gc(z)
GmT = mT(minv(G))
c1 = 1 / (h * (th[0] - th[1]) * (th[0] - th[2] + h) * (th[0] - th[3]))
print('c1 =', c1, c1 == Fr(603680, 25947))
def act_site(vecs, s, g):
    vs = list(vecs); v = vs[s]
    vs[s] = [sum(v[k] * g[k][j] for k in range(4)) for j in range(4)]
    return vs
for name, sites in [('vac', [unit(0), unit(3), unit(0), unit(3)]),
                    ('random product', [[Fr(1), Fr(2, 3), Fr(-1, 5), Fr(3)], [Fr(2), Fr(-1), Fr(1, 7), Fr(1, 2)], [Fr(-3), Fr(1), Fr(4, 9), Fr(1)], [Fr(1, 3), Fr(5), Fr(-2), Fr(1)]])]:
    Om = product_vec(sites)
    rhs = c1 * product_vec(act_site(sites, 0, GmT)).dot(T1)
    print('(i)', name, all(Om == rhs))

# (ii) twisted vacuum eigenbra
w = (Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9))
Vw = vandermonde(w)
vw = [w[0] ** 3, w[0] ** 2, w[0], Fr(1)]
uw = minv(Vw)[3]
chw = Chain(z=w)
for u in [Fr(3, 10), Fr(-7, 13)]:
    bra = product_vec([vw, uw, vw, uw])
    lhs = bra.dot(chw.t(1, u))
    ev = w[0] * (u - th[0] + h) * (u - th[1]) * (u - th[2] + h) * (u - th[3]) + (w[1] + w[2]) * math.prod(u - t for t in th) \
         + w[3] * (u - th[0]) * (u - th[1] - h) * (u - th[2]) * (u - th[3] - h)
    print('(ii) eigenbra at u=%s:' % u, all(lhs == ev * bra))
# (iii)
Gz = Gc(z)
q = sum(uw[i] * Gz[i][j] * vw[j] for i in range(4) for j in range(4))
pred = -math.prod(w[0] - zi for zi in z) / math.prod(w[3] - w[m] for m in range(3))
print('(iii)', q, pred, q == pred, q == Fr(232551, 1558000))

# frame map check t^{Gc} R = R t^{Lam}
V = vandermonde(z)
VmT = mT(minv(V))
def kron4(ms):
    out = np.array([[Fr(1)]], dtype=object)
    for m in ms:
        out = np.kron(out, np.array(m, dtype=object))
    return out
R = kron4([VmT, V, VmT, V])
chd = Chain(G=Gdiag(z))
for a, u in [(1, Fr(3, 10)), (2, Fr(-1, 3)), (3, Fr(5, 7))]:
    ok = (ch.t(a, u).dot(R) == R.dot(chd.t(a, u))).all()
    print('frame map a=%d' % a, ok)

# word basis exact
t0 = time.time()
Wd = sov_words(ch)
print('words built', time.time() - t0)
t0 = time.time()
print('rank mod p of 256 words:', rank_mod_p(Wd), time.time() - t0)
pickle.dump({'words': Wd, 'R': R}, open('words_exact.pkl', 'wb'))
