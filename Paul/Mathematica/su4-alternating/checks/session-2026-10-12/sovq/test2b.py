from sov_l4 import *
from fractions import Fraction as Fr
h, th, z = H0, TH0, Z0
ch = Chain(); T1 = ch.t(1, th[0])
GmT = mT(minv(Gc(z)))
c1 = 1 / (h * (th[0] - th[1]) * (th[0] - th[2] + h) * (th[0] - th[3]))
v = [Fr(1), Fr(2, 3), Fr(-1, 5), Fr(3)]
r = [Fr(2), Fr(-1), Fr(1, 7), Fr(1, 2)]
vr = sum(a*b for a, b in zip(v, r)); vv = sum(a*a for a in v)
w = [ri - vr/vv*vi for ri, vi in zip(r, v)]
sites = [v, w, v, w]
Om = product_vec(sites)
s1 = [sum(v[k]*GmT[k][j] for k in range(4)) for j in range(4)]
rhs = c1 * product_vec([s1, w, v, w]).dot(T1)
print('(i) BPS v.w=0:', all(Om == rhs))
