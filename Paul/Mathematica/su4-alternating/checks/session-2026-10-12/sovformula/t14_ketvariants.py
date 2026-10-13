"""Find a full-rank regularised KET SoV basis at the ABJM point (L = 2), variants of the mirror of the bra regularisation."""
import time
from regbases import *
R = Reg()
h = R.h; th = R.th
ops = []
for p in (0, 1):
    thA = th[2 * p]; uF = th[2 * p + 1]; u0 = uF + h
    t1A, d1A = R.dop(1, thA); t1u0, d1u0 = R.dop(1, u0); t2A, d2A = R.dop(2, thA); t3, d3 = R.dop(3, thA + h)
    ops.append(dict(t1A=t1A, t1u0=t1u0, Wa=d2A - d1A.dot(t1u0), Wb=d2A - t1A.dot(d1u0), Wc=d2A - d1A.dot(t1u0) - t1A.dot(d1u0), t2A=t2A, Z=d3))
vac = S.product_vec([S.unit(i) for i in S.KET_VAC])
def kets(wname, special3):
    cols = []
    for (m1, n2, m3, n4) in S.WORDS:
        v = vac
        for p, (m, n) in reversed(list(enumerate(((m1, n2), (m3, n4))))):
            P = ops[p]
            if n == 1: v = P['t1u0'].dot(v)
            elif n == 2: v = (P['t2A'] if (special3 and m == 3) else P[wname]).dot(v)
            elif n == 3: v = P['Z'].dot(v)
            for _ in range(m): v = P['t1A'].dot(v)
        cols.append(v)
    return np.array(cols).T
for wname in ('Wa', 'Wb', 'Wc'):
    for sp3 in (False, True):
        t0 = time.time(); Y = kets(wname, sp3)
        print(wname, 'special m=3 ->', sp3, ': rank', S.rank_mod_p(Y.T), ' (%.0f s)' % (time.time() - t0), flush=True)
