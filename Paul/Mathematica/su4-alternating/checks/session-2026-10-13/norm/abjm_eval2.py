"""Python fallback: eigenvalue data, 1/N = sum tauT M tau at the ABJM point (eps-chain), direct check, and Q-based C_Q.
usage: python abjm_eval.py 1/100 100   (needs gram_1_100.pkl from gram_abjm.py)"""
import sys, time, pickle, json
SOV = r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/checks/session-2026-10-12/sovformula'
sys.path.insert(0, SOV)
import numpy as np, mpmath as mp
from fractions import Fraction as Fr
import sov_l4 as S
from mpu import Ctx, tomp
eps = Fr(sys.argv[1]); dps = int(sys.argv[2]); mp.mp.dps = dps
tag = str(eps).replace('/', '_')
d = pickle.load(open('gram_%s.pkl' % tag, 'rb'))
G, ops, z, thF = d['G'], d['ops'], d['z'], d['thF']
phi = [Fr(9, 10), Fr(7, 20), Fr(3, 5), Fr(-11, 10)]; thdir = [Fr(0), Fr(37, 100)]; c = Fr(3, 10); h = Fr(1, 2)
th = (thF[0] + 2 * h, thF[0], thF[1] + 2 * h, thF[1])
t0 = time.time()
ctx = Ctx(2, eps, phi, thdir, dps=dps)
w = (0, 0, 0, 0); sec = ctx.sectors[w]; ns = len(sec)
A1 = ctx.secmat(w, 1, Fr(3, 11)); A3 = ctx.secmat(w, 3, Fr(-5, 7)); A2 = ctx.secmat(w, 2, Fr(2, 9))
Mg = [[A1[i][j] + Fr(7, 10) * A3[i][j] - Fr(3, 13) * A2[i][j] for j in range(ns)] for i in range(ns)]
Am = mp.matrix([[tomp(x) for x in r] for r in Mg])
E, ER = mp.eig(Am); EL = mp.inverse(ER)
print('eig %.0f s' % (time.time() - t0), flush=True)
V = S.vandermonde(list(z)); Vi = S.minv(V); VmT = S.mT(Vi); VT = S.mT(V)
def nzrow(op, p): return [(j, tomp(op[p, j])) for j in range(256) if op[p, j] != 0]
# Gram -> mp and LU (cached by mpmath)
Gm = mp.matrix(256, 256)
for i in range(256):
    for j in range(256):
        if G[i, j] != 0: Gm[i, j] = tomp(G[i, j])
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
lam = [Fr(1)/ (lambda zz, j: __import__('math').prod([zz[j]-zz[m] for m in range(4) if m != j]))(z, j) for j in range(4)]
aa = (phi[0] - c) / (phi[0] - phi[3]); kap = 1 - c * eps; y = (z[1] + z[2]) / 2; s = 1 - kap * (z[0] + z[3]); Nk = (1 - kap * z[0]) ** 2
kap0 = -((z[0] - z[1]) * (z[0] - z[2])) / ((z[3] - z[1]) * (z[3] - z[2])); nrm = lam[0] ** 2 * ((z[3] - z[1]) * (z[3] - z[2])) ** 2
nu = th[1] - th[3]; ut = [1, -(z[1] + z[2]), z[1] * z[2]]; u0 = th[1] + h
def lamK(k, u):
    f = Fr(1)
    for x in [u - th[2] - j * h for j in range(k)] + [u - th[3] - k * h] + [u - th[3] - j * h for j in range(k - 1)]:
        if x != 0: f *= x
    return f
l1, l2, l3 = lamK(1, th[1]), lamK(2, u0), lamK(3, u0)
def Pab(r2, r4):
    nu_, h_ = tomp(nu), tomp(h)
    return [[(nu_/(nu_+h_))**2*r2**2*r4**2, nu_/(nu_+h_)*r2**2*r4, (nu_-h_)/nu_*r2**2 + h_/(nu_+h_)*r2*r4],
            [nu_/(nu_+h_)*r2*r4**2 - h_/(nu_+h_)*r2**2*r4, nu_/(nu_+h_)*r2*r4, r2],
            [r4**2 - (h_/nu_)*(1+h_/(nu_+h_))*r2*r4 - (h_/nu_)*(nu_-h_)/nu_*r2**2, r4 - (h_/nu_)*r2, 1]]
bra0W = ctx.prodbra([[1, 0, 0, aa], [-aa, 0, 0, 1], [1, 0, 0, aa], [-aa, 0, 0, 1]], w)
bra0D = ctx.prodbra([lam, [1]*4, lam, [1]*4], w); ket0D = ctx.prodbra([[1]*4, lam, [1]*4, lam], w)
rows = []
LUdone = False
for m in range(ns):
    psi = [ER[i, m] for i in range(ns)]; psiL = [EL[m, i] for i in range(ns)]
    nv = mp.sqrt(mp.fsum(abs(x) ** 2 for x in psi))
    vd = [mp.mpf(0)] * 256; ld = [mp.mpf(0)] * 256
    for r, cidx in enumerate(sec): vd[cidx] = psi[r]; ld[cidx] = psiL[r]
    Pc = S._apply_kron([VmT, V, VmT, V], vd); Lc = S._apply_kron([VT, Vi, VT, Vi], ld, left=True)
    p = max(range(256), key=lambda i: abs(Pc[i]))
    def ev(op):
        return mp.fsum(x * Pc[j] for j, x in nzrow(op, p)) / Pc[p]
    br, kt = ops['bra'], ops['ket']
    F = [ev(br[q]['t1F']) for q in (0, 1)]
    Ob = [[1] + [ev(br[q]['O'][a]) for a in (1, 2, 3)] for q in (0, 1)]; O23 = [ev(br[q]['t2u0']) for q in (0, 1)]
    Av = [ev(kt[q]['t1A']) for q in (0, 1)]; Kv = [[1] + [ev(kt[q]['K'][n]) for n in (1, 2, 3)] for q in (0, 1)]
    tauX = []; tauY = []
    for (a1, k2, a3, k4) in S.WORDS:
        f = []
        for q, (a, k) in enumerate(((a1, k2), (a3, k4))):
            o = (O23[q] if (a == 2 and k == 3) else Ob[q][a]) * F[q] ** k
            f.append(o)
        tauX.append(f[0] * f[1])
        g = [Av[q] ** mm * Kv[q][nn] for q, (mm, nn) in enumerate(((a1, k2), (a3, k4)))]
        tauY.append(g[0] * g[1])
    sol = mp.lu_solve(Gm, mp.matrix(tauX))
    invNQ = mp.fsum(tauY[i] * sol[i] for i in range(256))
    ov0 = Pc[51]; ovL0 = Lc[204]; pr = mp.fsum(Lc[i] * Pc[i] for i in range(256)); invNd = pr / (ov0 * ovL0)
    pu = [x / nv for x in psi]
    realQ = abs(E[m].imag) < mp.mpf(10) ** (-dps // 2)
    ov0D = mp.fsum(a * b for a, b in zip(bra0D, pu)); wDir = mp.fsum(a * b for a, b in zip(bra0W, pu))
    X1, r4, X2, X3 = F[0], F[1], O23[0], Ob[0][3]
    K = s**2 + tomp(kap)*tomp(s)*(X1/(tomp(h)**2*tomp(l1)) + X2/(2*tomp(h)**4*tomp(y)*tomp(l2)) + X3/(6*tomp(h)**5*tomp(y)**2*tomp(l3))) + tomp(kap)**2*X3*X1/(24*tomp(h)**7*tomp(y)**2*tomp(l1)*tomp(l3))
    rho2 = X1 / (tomp(h) * tomp(th[0] - th[1]) * tomp(th[1] - th[2]) * tomp(nu - h)); rho4 = r4 / (tomp(h) * tomp(th[0] - th[3]) * tomp(nu) * tomp(th[2] - th[3]))
    P = Pab(rho2, rho4); utm = [tomp(x) for x in ut]
    utPut = mp.fsum(utm[i] * P[i][j] * utm[j] for i in range(3) for j in range(3))
    fac = (tomp(aa) / tomp(kap0)) ** 2 * utPut / tomp(nrm) / tomp(Nk)
    Cdir = abs(K * wDir / tomp(Nk)); CQ = abs(K * fac) * mp.sqrt(abs(1 / invNQ)); Hn = abs(ov0D) ** 2 / (abs(1 / invNd) * mp.fsum(abs(x) ** 2 for x in pu))
    rows.append(dict(m=m, real=bool(realQ), Cdir=Cdir, CQ=CQ, invNerr=abs(invNQ / invNd - 1), H=Hn, data=dict(F=F, Ob=Ob, O23=O23, Av=Av, Kv=Kv, E=E[m])))
    print(m, realQ, mp.nstr(Cdir, 12), mp.nstr(CQ, 12), mp.nstr(abs(invNQ/invNd - 1), 3), mp.nstr(Hn, 10), '%.0f s' % (time.time() - t0), flush=True)
pickle.dump(rows, open('evaldata_%s.pkl' % tag, 'wb'))
