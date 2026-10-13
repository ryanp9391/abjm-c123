"""High-precision SoV side of the untwisting check: sum_x c_x(b) W_x <0|Psi> vs direct <B|Psi> vs data C2.
Exact c (multi-prime CRT), exact weight-sector blocks of the diagonal-frame transfer matrices (multi-prime CRT), mpmath eigenvectors."""
import numpy as np, sys, time, pickle, math
import mpmath as mp
from sympy import prevprime
from core import *
from exact import crt_pair, ratrec_big
from nb_common import WT
t0 = time.time()
eps = Fr(1, int(sys.argv[1])); dps = int(sys.argv[2]); NPR = int(sys.argv[3])
mp.mp.dps = dps
bval = Fr(7, 10)
phi = [Fr(9, 10), Fr(2, 5), Fr(-3, 10), Fr(-11, 10)]; thA = [Fr(1), Fr(0), Fr(1), Fr(0)]
th = [thA[s] + eps * th0[s] for s in range(4)]; z = [1 + eps * phi[j] for j in range(4)]
h = h0
primes = []; q = 2 ** 26
while len(primes) < NPR:
    q = prevprime(q); primes.append(q)
secs = {}
for k, w in enumerate(WT): secs.setdefault(w, []).append(k)
lam = [Fr(1) / math.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
import sympy as sp
Vs = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in Vand(z)]); Vis = Vs.inv()
vz4 = [z[3] ** 3, z[3] ** 2, z[3], Fr(1)]
l1 = [Fr(int(sp.fraction(Vis[0, j])[0]), int(sp.fraction(Vis[0, j])[1])) for j in range(4)]
Bc_sites = [[Fr(1) if k == 0 else Fr(0) for k in range(4)], None, [1, 0, 0, 0], [0, 0, 0, 1]]
A_vec = [Fr(int(k == 0)) + bval * lam[0] * vz4[k] for k in range(4)]
F_vec = [Fr(int(k == 3)) - bval * l1[k] for k in range(4)]
Gd = [[z[i] if i == j else 0 for j in range(4)] for i in range(4)]
ops = [('H1', 1, Fr(37, 100)), ('H2', 2, Fr(-113, 100)), ('H3', 3, Fr(3, 2)), ('H4', 1, Fr(-113, 100)),
       ('A1_1', 1, th[0] - h), ('A1_2', 2, th[0] - h), ('A1_3', 3, th[0] - h), ('F2', 1, th[1]),
       ('A3_1', 1, th[2] - h), ('A3_2', 2, th[2] - h), ('A3_3', 3, th[2] - h), ('F4', 1, th[3])]
blocks_mod = {name: {} for name, _, _ in ops}
cmod = {}
for pp in primes:
    R = Ring('mod', pp)
    chd = Chain(R, th=th, G=Gd)
    for name, a, u in ops:
        T = chd.t(a, u)
        blocks_mod[name][pp] = {w: T[np.ix_(idx, idx)].copy() for w, idx in secs.items()}
    chc = Chain(R, th=th)  # companion, z0 (words are z-independent)
    M, vac = sov_words(chc)
    Y = R.inv(M)
    bvec = kronvecs(R, [R.arr(A_vec), R.arr(F_vec), R.arr([1, 0, 0, 0]), R.arr([0, 0, 0, 1])])
    cmod[pp] = R.mm(bvec[None, :], Y)[0]
    print('prime done', pp, round(time.time() - t0, 1), flush=True)
def rec(vals):
    a_, m_ = int(vals[0]) % primes[0], primes[0]
    for v, pp in zip(vals[1:-1], primes[1:-1]): a_, m_ = crt_pair(a_, m_, int(v) % pp, pp)
    r = ratrec_big(a_, m_)
    if r is None or (r.numerator - r.denominator * int(vals[-1])) % primes[-1] != 0: return None
    return r
# exact c
c = {}; fails = 0
for x in range(256):
    vals = [cmod[pp][x] for pp in primes]
    if all(v == 0 for v in vals[:3]): continue
    r = rec(vals)
    if r is None: fails += 1
    c[x] = r
print('c: nonzero', len(c), 'reconstruction failures', fails, flush=True)
if fails: sys.exit('increase primes')
# exact blocks
blocks = {name: {} for name, _, _ in ops}; bf = 0
for name, _, _ in ops:
    for w, idx in secs.items():
        n = len(idx); B = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                vals = [blocks_mod[name][pp][w][i, j] for pp in primes]
                if all(v == 0 for v in vals[:3]): B[i][j] = Fr(0); continue
                r = rec(vals)
                if r is None: bf += 1; r = Fr(0)
                B[i][j] = r
        blocks[name][w] = B
print('t-blocks reconstructed, failures', bf, round(time.time() - t0, 1), flush=True)
if bf: sys.exit('increase primes')
tomp = lambda B: mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in B])
# data
import wxf
D = wxf.load(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf')
def sp2d(sa):
    dims = int(sa.args[1][0]); cols = sa.args[3][1][1]; vals = sa.args[3][2]
    v = np.zeros(dims, dtype=complex)
    for cc, xx in zip(np.array(cols).ravel(), vals): v[int(cc) - 1] = xx
    return v
data = [(tuple(int(x) for x in s['w']), [int(x) for x in s['dyn']], s['E'], s['C'], sp2d(s['psi'])) for s in D['L2']]
del D
lamm = [mp.mpf(x.numerator) / x.denominator for x in lam]
Bdiag_sites = [[lamm[k] + (mp.mpf(bval.numerator) / bval.denominator) * lamm[0] * (k == 3) for k in range(4)],
               [1 - (mp.mpf(bval.numerator) / bval.denominator) * (k == 0) for k in range(4)], lamm, [mp.mpf(1)] * 4]
vac0_sites = [lamm, [mp.mpf(1)] * 4, lamm, [mp.mpf(1)] * 4]
def site_comp(sites, idx):
    i = wordsIdx[idx]; return sites[0][i[0]] * sites[1][i[1]] * sites[2][i[2]] * sites[3][i[3]]
coefH = [mp.mpf(1), mp.mpf('-0.7'), mp.mpf('0.3'), mp.mpf('1.1')]
print('| ε = δ | Dynkin | E_MZ | # c_x ≠ 0 | SoV Σc_xW_x⟨0|Ψ⟩/(λ₁²(−b)^{2−w₁}) | direct ⟨B|Ψ⟩/(…) | rel. diff | data C₂ | \|⟨0|Ψ⟩\|/λ₁² |')
for (w, dyn, E, C, psi) in data:
    idx = secs[w]
    Hb = sum((coefH[k] * tomp(blocks[name][w]) for k, name in enumerate(['H1', 'H2', 'H3', 'H4'])), mp.zeros(len(idx), len(idx)))
    ev, U = mp.eig(Hb)
    # match by overlap with data psi (float)
    best, bov = None, -1
    for k in range(len(idx)):
        vk = np.array([complex(U[i, k]) for i in range(len(idx))])
        ov = abs(np.vdot(psi[idx], vk)) / np.linalg.norm(vk)
        if ov > bov: best, bov = k, ov
    v = [U[i, best] for i in range(len(idx))]
    nrm = mp.sqrt(sum(abs(x) ** 2 for x in v)); v = [x / nrm for x in v]
    ph = sum(mp.mpc(complex(psi[idx[i]]).real, -complex(psi[idx[i]]).imag) * v[i] for i in range(len(idx)))
    v = [x * abs(ph) / ph for x in v]
    def tau(name):
        Bm = tomp(blocks[name][w]); tv = Bm * mp.matrix(v)
        return sum(mp.conj(v[i]) * tv[i] for i in range(len(idx))) / sum(mp.conj(v[i]) * v[i] for i in range(len(idx)))
    A1 = [1, tau('A1_1'), tau('A1_2'), tau('A1_3')]; A3 = [1, tau('A3_1'), tau('A3_2'), tau('A3_3')]
    r2 = tau('F2'); r4 = tau('F4')
    K = sum((mp.mpf(cv.numerator) / cv.denominator) * A1[wordsIdx[x][0]] * r2 ** wordsIdx[x][1] * A3[wordsIdx[x][2]] * r4 ** wordsIdx[x][3] for x, cv in c.items())
    ov0 = sum(site_comp(vac0_sites, idx[i]) * v[i] for i in range(len(idx)))
    direct = sum(site_comp(Bdiag_sites, idx[i]) * v[i] for i in range(len(idx)))
    nrmz = lamm[0] ** 2 * (-(mp.mpf(bval.numerator) / bval.denominator)) ** (2 - w[0])
    sov = K * ov0 / nrmz; dr = direct / nrmz
    rd = abs(sov - dr) / max(abs(dr), mp.mpf(10) ** (-dps + 5))
    f = lambda x: mp.nstr(mp.re(x), 10) + ('' if abs(mp.im(x)) < mp.mpf(10) ** -30 else ' ' + mp.nstr(mp.im(x), 6) + 'i')
    print(f'| {float(eps):.0e} | {dyn} | {E:.4f} | {len(c)} | {f(sov)} | {f(dr)} | {mp.nstr(rd, 3)} | {complex(C[2]).real:.10f} | {mp.nstr(abs(ov0) / lamm[0] ** 2, 3)} |', flush=True)
print('time', round(time.time() - t0, 1))
