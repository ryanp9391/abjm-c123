import numpy as np
from core import *

R = Ring('mod')
ch = Chain(R)
h, th = ch.h, ch.th
t1i = R.inv(ch.t(1, th[0]))
GmT = R.inv(R.arr(ch.G)).T.copy()
c1 = Fr(1) / (h * (th[0] - th[1]) * (th[0] - th[2] + h) * (th[0] - th[3]))


def test(vs, name):
    Om = kronvecs(R, vs)
    lhs = R.mm(Om[None, :], t1i)[0]
    v1r = R.mm(R.arr(vs[0])[None, :], GmT)[0]
    rhs = kronvecs(R, [v1r] + vs[1:])
    res = R.red(lhs - R.smul(c1, rhs))
    print(name, 'residual nonzero count', np.count_nonzero(res))


v = [1, 2, -1, 3]; u = [2, -1, 1, 1]  # v.u = 2-2-1+3 = 2
u = [2, -1, 1, -1]  # v.u = 2-2-1-3=-4
u = [1, 1, 0, -1]   # 1+2+0-3 = 0
test([v, u, v, u], 'global BPS')
v2 = [3, 1, 1, 2]; u2 = [1, -1, 0, -1]  # 3-1+0-2=0
test([v, u, v2, u2], 'two-block BPS')
test([v, [1, 2, 3, 4], v2, u2], 'pair1 non-BPS')
test([v, u, v2, [1, 2, 3, 4]], 'pair2 non-BPS')
test([v, [1, 2, 3, 4], v2, [2, 1, 1, 1]], 'both non-BPS')
