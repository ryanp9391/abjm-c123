"""obstruction.py -- (a) charged-sector eigenvalues of the five words; (b) generic-z (z2 != z3) consistent test of the formal twisted vacuum;
(c) best-fit residual in the untwisting scaling.  Python fallback for notebook sections 5d, 5e.  Output: obstruction_out.txt"""
from rm_base import *
from helpers import tw_form, tw_u, prod_, pred, comp_v
import io, contextlib
out = open('obstruction_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

# (a) charged eigenvalues, exact, two random generic z
for zz in ([Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)], [Fr(-1, 3), Fr(7, 5), Fr(2, 9), Fr(11, 4)]):
    ops = pair_ops_z(zz); sig = sum(zz); det = prod_(zz)
    norm = [1 / (2 * h**2), -1 / (2 * h**4), 1 / (12 * h**5)]
    allok = True
    for a in range(4):
        for b in range(4):
            if a == b: continue
            r = [Fr(0)] * 16; r[a * 4 + b] = Fr(1)
            pred = [zz[b], zz[b] * (sig - zz[a] - zz[b]), det / zz[a]]
            for k in range(3):
                s = [x * norm[k] for x in matvec(r, ops[k])]
                eig = s[a * 4 + b]
                nz = [i for i in range(16) if s[i] != 0 and i != a * 4 + b]
                allok &= (not nz) and eig == pred[k]
    P('(a) charged weight states e_a (x) e_b (a != b) are eigenvectors of x1, x2, x3 with eigenvalues z_b, z_b(sigma - z_a - z_b), det/z_a:', allok, ' z =', [str(x) for x in zz])

# (b) consistent generic-z test of the formal twisted vacuum, and of W_a
zg = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]; ww = [Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)]
def obs(name, v, u, zz):
    ops = pair_ops_z(zz); pw = words(v, u, ops); comps = kap_comps(v, u)
    S = [(a, b) for a in range(4) for b in range(4) if a != b and v[a] * u[b] != 0]
    sig = sum(zz); det = prod_(zz)
    M = [[Fr(1), zz[b], zz[b] * (sig - zz[a] - zz[b]), det / zz[a], det * zz[b] / zz[a]] for (a, b) in S]
    rk = rank(pw); rc = rank(comps); rj = rank(pw + comps)
    kaps = KAPS + [1 / zz[3]]
    rkT = [rank(pw + [target(v, u, k, zz)]) for k in kaps]
    P(name, '| v.u =', sum(x * y for x, y in zip(v, u)), '| rank span', rk, '| #kappa-comps', len(comps), 'rank', rc, '| dim(span cap kappa-orbit span) =', rk + rc - rj,
      '| charged |S| =', len(S), 'rank of |S|x5 eigenvalue matrix', rank(M) if M else 0, '| rank(span+T_kappa)', rkT)
vf, uf = tw_form(ww, zg)
obs('(b) formal twisted vacuum, generic z=(3/2,5/4,2/3,4/7)', vf, uf, zg)
a = Fr(1, 3)
obs('(b) W_a, a=1/3, same generic z', vec(1, 0, 0, a), vec(-a, 0, 0, 1), zg)

# (c) scaling residuals
P('(c) best-fit relative residual of T_kappa in the 5-word span; z = 1 + eps*phi, w = 1 + eps*omega, kappa = 1 - c*eps; phi = (9/10, 7/20, 7/20 + dphi, -11/10), omega = (7/10,-3/10,1/5,1/2), c = 4/5')
charged = [a * 4 + b for a in range(4) for b in range(4) if a != b]
P('dphi  eps | tw.vac all | tw.vac charged | W_a all | W_a charged')
for dl in (Fr(1, 100), Fr(1, 10), Fr(1, 2)):
    for e in (2, 3, 4):
        ep = Fr(1, 10**e); phi = [Fr(9, 10), Fr(7, 20), Fr(7, 20) + dl, Fr(-11, 10)]; om = [Fr(7, 10), Fr(-3, 10), Fr(1, 5), Fr(1, 2)]
        zz = [1 + ep * p for p in phi]; w_ = [1 + ep * o for o in om]; kap = 1 - Fr(4, 5) * ep
        v, u = tw_form(w_, zz); ops = pair_ops_z(zz)
        pw = words(v, u, ops); T = target(v, u, kap, zz)
        wv, wu = vec(1, 0, 0, Fr(1, 3)), vec(Fr(-1, 3), 0, 0, 1)
        pw2 = words(wv, wu, ops); T2 = target(wv, wu, kap, zz)
        P(float(dl), 'e=%d' % e, '|', '%.4g' % residual_ls(pw, T), '%.4g' % residual_ls(pw, T, charged), '|', '%.4g' % residual_ls(pw2, T2), '%.4g' % residual_ls(pw2, T2, charged))
out.close()
