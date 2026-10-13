"""Try right (ket) SoV bases |y> = word(t)|Omega>; measure N = <x|y>; sparsity."""
import numpy as np, itertools
from core import *
R = Ring('mod', P1); p = P1; ch = Chain(R); h = ch.h; th = ch.th
M, vac = sov_words(ch)
I = np.eye(256, dtype=np.int64)
def mpow(A, k):
    out = I
    for _ in range(k): out = R.mm(out, A)
    return out
def ket_words(omega, Aop, Fop):
    """|y> = A-site ops on sites 1,3 and F-site ops on sites 2,4 acting on ket omega: columns"""
    cols = []
    A1, A3 = Aop(0), Aop(2); F2, F4 = Fop(1), Fop(3)
    for (a1, k2, a3, k4) in wordsIdx:
        v = omega[:, None]
        v = R.mm(A1[a1], v); v = R.mm(F2[k2], v); v = R.mm(A3[a3], v); v = R.mm(F4[k4], v)
        cols.append(v[:, 0])
    return np.array(cols).T
variants = {}
for om_name, om in [('e4 e1 e4 e1', kronvecs(R, [e(4), e(1), e(4), e(1)])), ('e1 e4 e1 e4', vac)]:
    for apt in (0, h, -h):
        for fsh in (h, 0, -h, 2 * h):
            Aop = lambda s, apt=apt: [mpow(ch.t(1, th[s] + apt), k) for k in range(4)]
            Fop = lambda s, fsh=fsh: [ch.t(a, th[s] + fsh) for a in range(4)]
            Yk = ket_words(om, Aop, Fop)
            r = modrank(Yk, p)
            if r < 256:
                print(om_name, 'A: t1(thA%+.1f)^k' % float(apt), 'F: t_a(thF%+.1f)' % float(fsh), 'rank', r); continue
            N = R.mm(M, Yk)
            print(om_name, 'A: t1(thA%+.1f)^k' % float(apt), 'F: t_a(thF%+.1f)' % float(fsh), 'rank 256, nnz <x|y> =', np.count_nonzero(N))
