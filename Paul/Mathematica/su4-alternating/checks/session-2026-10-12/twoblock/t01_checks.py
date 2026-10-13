import time, numpy as np
from core import *

t0 = time.time()
for mode in ['float', 'mod']:
    R = Ring(mode)
    ch = Chain(R)
    h, th = ch.h, ch.th
    u1, u2 = Fr(3, 10), Fr(-7, 10)
    A = ch.t(1, u1); B3 = ch.t(3, u2); B2 = ch.t(2, u2)
    c13 = R.red(R.mm(A, B3) - R.mm(B3, A)); c12 = R.red(R.mm(A, B2) - R.mm(B2, A))
    c23 = R.red(R.mm(ch.t(2, u1), B3) - R.mm(B3, ch.t(2, u1)))
    print(mode, 'commutators max', np.abs(c13).max(), np.abs(c12).max(), np.abs(c23).max(), 'time', time.time() - t0)
    # calibration (i): <Omega| t1(th1)^{-1} = c1 <Omega| rho1(G^{-T})
    vs = [[1, 2, -1, 3], [2, -1, 1, 1], [-1, 1, 2, 3], [1, 1, -2, 1]]
    Om = kronvecs(R, vs)
    t1i = R.inv(ch.t(1, th[0]))
    lhs = R.mm(Om[None, :], t1i)[0]
    GmT = R.inv(R.arr(ch.G)).T.copy()
    v1r = R.mm(R.arr(vs[0])[None, :], GmT)[0]
    rhs = kronvecs(R, [v1r] + vs[1:])
    c1 = Fr(1) / (h * (th[0] - th[1]) * (th[0] - th[2] + h) * (th[0] - th[3]))
    print('c1 =', c1)
    print(mode, 'calib (i) residual', np.abs(R.red(lhs - R.smul(c1, rhs))).max())
    # calibration (ii): twisted-vacuum eigenbra of t1^{Gc(w)}
    w = [Fr(7, 3), Fr(-2, 5), Fr(5, 6), Fr(9, 7)]
    chw = Chain(R, z=w)
    vw = [w[0] ** 3, w[0] ** 2, w[0], 1]
    Vw = np.array(Vand(w), dtype=object)
    import sympy
    Vi = sympy.Matrix(Vw.tolist()).inv()
    uw = [Fr(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1])) for x in Vi.row(3)]
    Bw = kronvecs(R, [vw, uw, vw, uw])
    u = Fr(13, 17)
    ev = w[0] * (u - th[0] + h) * (u - th[1]) * (u - th[2] + h) * (u - th[3]) + (w[1] + w[2]) * np.prod([u - x for x in th]) + w[3] * (u - th[0]) * (u - th[1] - h) * (u - th[2]) * (u - th[3] - h)
    res = R.red(R.mm(Bw[None, :], chw.t(1, u))[0] - R.smul(ev, Bw))
    print(mode, 'calib (ii) residual', np.abs(res).max())
    M, vac = sov_words(ch)
    print(mode, 'rank words', R.rank(M), 'time', time.time() - t0)
