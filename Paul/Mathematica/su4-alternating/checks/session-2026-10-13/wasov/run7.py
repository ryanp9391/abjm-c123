"""Sectors, exactness of the component SoV sums on 256 eigenstates (mp 50), eigenbra property of the n = 0, 4 components."""
from wa_sov_core import *
import numpy as np
import sov_l4 as S
import mpmath as mp
mp.mp.dps = 50
ch, W, Y, Yk, G, M = build(verbose=False)
eig = cached('eigen', S.eigenstates)
sts = [S.Eig(d) for d in eig]
right, left = cached('lr', lambda: S.companion_vectors(eig))
cn = pickle.load(open(os.path.join(CACHE, 'cn.pkl'), 'rb'))
V = vandermonde(Z0); Vi = minv(V)
vz = [[V[i][j] for i in range(4)] for j in range(4)]; ut = [Vi[j] for j in range(4)]
comp = {n: np.array([Fr(0)] * 256, dtype=object) for n in range(5)}
for flips in itertools.product([0, 1], repeat=4):
    vs = [((vz[3] if fl else vz[0]) if s % 2 == 0 else ([-x for x in ut[0]] if fl else ut[3])) for s, fl in enumerate(flips)]
    comp[sum(flips)] = comp[sum(flips)] + product_vec(vs)
def sect(w):
    return 2 - w[0] if (w[1] == 0 and w[2] == 0 and w[0] + w[3] == 0 and abs(w[0]) <= 2) else None
print('sector sizes:', {n: sum(1 for st in sts if sect(st.weight) == n) for n in range(5)}, ' states outside:', sum(1 for st in sts if sect(st.weight) is None))
idx = {w: i for i, w in enumerate(WORDS)}
err = {n: mp.mpf(0) for n in range(5)}; inmax = {n: mp.mpf(0) for n in range(5)}; offmax = {n: mp.mpf(0) for n in range(5)}
G1 = {}
for k, st in enumerate(sts):
    tau_site = [[(mp.mpf(1) if lab == 0 else (st.tau(lab, S._th[s] - S._h) if S.TYPES[s] == 'A' else st.tau(1, S._th[s]) ** lab)) for lab in range(4)] for s in range(4)]
    tx = [mp.fprod([tau_site[s][x[s]] for s in range(4)]) for x in WORDS]
    ov0 = right[k][51]
    for n in range(5):
        direct = mp.fsum(S.tomp(comp[n][c]) * right[k][c] for c in range(256) if comp[n][c] != 0)
        form = mp.fsum(S.tomp(cn[n][i]) * tx[i] for i in range(256)) * ov0
        err[n] = max(err[n], abs(direct - form))
        if sect(st.weight) == n:
            inmax[n] = max(inmax[n], abs(direct)); G1[(n, k)] = direct / ov0
        else:
            offmax[n] = max(offmax[n], abs(direct))
for n in range(5):
    print('component', n, ' max|direct - SoV sum| =', mp.nstr(err[n], 3), ' max in-sector |overlap| =', mp.nstr(inmax[n], 5), ' max off-sector |overlap| =', mp.nstr(offmax[n], 3))
# eigenbra property of components 0 and 4
chz = Chain()
for n, perm in ((0, (0, 1, 2, 3)), (4, (3, 1, 2, 0))):
    for u in (Fr(3, 10), Fr(-7, 10)):
        op = chz.t(1, u); b = comp[n]
        r = b.dot(op); p = next(i for i, x in enumerate(b) if x != 0); e = r[p] / b[p]
        prop = all(r[i] == e * b[i] for i in range(256))
        zz = [Z0[i] for i in perm]
        Vv = lambda i, x: S.Vvac(i, x, h=H0, th=TH0)
        ev = sum(zz[i] * S.V4_exact(u) if False else 0 for i in range(0))
        print('component', n, 'u =', u, 'proportional to eigenbra of t1:', prop, ' eigenvalue', float(e))
