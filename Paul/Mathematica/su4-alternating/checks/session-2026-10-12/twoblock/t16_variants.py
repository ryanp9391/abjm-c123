"""Word-basis variants: F-site telescoping products, A-site variants; sparsity of physical bra & families."""
import numpy as np, sympy
from core import *
R = Ring('mod', P1); p = P1; ch = Chain(R); h = ch.h; th = ch.th
z = z0
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
V = R.arr(Vand(z)); Vi = modinv(V, p); VT = V.T.copy()
def vm(v, A): return R.mm(R.arr(v)[None, :] if not isinstance(v, np.ndarray) else v[None, :], A)[0]
e1, e4 = R.arr(e(1)), R.arr(e(4))
vz4 = vm(e(4), VT); l1 = vm(e(1), Vi)
lam1 = R.c(lam[0])
def phys(b):
    b = R.c(b)
    return [(e1 + b * lam1 % p * vz4) % p, (e4 - b * l1) % p, e1, e4]
I = np.eye(256, dtype=np.int64)
def Fops(s, mode):
    t = lambda u: ch.t(1, u)
    if mode == 'pow': return [I, t(th[s]), R.mm(t(th[s]), t(th[s])), R.mm(R.mm(t(th[s]), t(th[s])), t(th[s]))]
    if mode == 'tele+': return [I, t(th[s]), R.mm(t(th[s]), t(th[s] + h)), R.mm(R.mm(t(th[s]), t(th[s] + h)), t(th[s] + 2 * h))]
    if mode == 'tele-': return [I, t(th[s]), R.mm(t(th[s]), t(th[s] - h)), R.mm(R.mm(t(th[s]), t(th[s] - h)), t(th[s] - 2 * h))]
def Aops(s, mode):
    if mode == 'ta': return [ch.t(a, th[s] - h) for a in range(4)]
    if mode == 't3pow':
        t3 = ch.t(3, th[s] - h); return [I, t3, R.mm(t3, t3), R.mm(R.mm(t3, t3), t3)]
    if mode == 't1pow':
        t1 = ch.t(1, th[s] - h); return [I, t1, R.mm(t1, t1), R.mm(R.mm(t1, t1), t1)]
    if mode == 'ta_tele':  # t_a fused at theta-h with telescoping t3
        t3 = ch.t(3, th[s] - h); return [I, t3, R.mm(t3, ch.t(3, th[s] - 2 * h)), R.mm(R.mm(t3, ch.t(3, th[s] - 2 * h)), ch.t(3, th[s] - 3 * h))]
def words(amode, fmode):
    vac = kronvecs(R, [e1, e4, e1, e4])
    A1, A3 = Aops(0, amode), Aops(2, amode); F2, F4 = Fops(1, fmode), Fops(3, fmode)
    rows = []
    for (a1, k2, a3, k4) in wordsIdx:
        v = R.mm(vac[None, :], A1[a1]); v = R.mm(v, F2[k2]); v = R.mm(v, A3[a3]); v = R.mm(v, F4[k4]); rows.append(v[0])
    return np.array(rows)
rng = np.random.default_rng(3)
def bps():
    v = rng.integers(1, p, 4); u = rng.integers(1, p, 4); s = int((v[:3] * u[:3] % p).sum() % p); u[3] = (-s * pow(int(v[3]), -1, p)) % p; return v, u
v, u = bps()
for am in ('ta', 't3pow', 't1pow', 'ta_tele'):
    for fm in ('pow', 'tele+', 'tele-'):
        M = words(am, fm); r = modrank(M, p)
        if r < 256: print(am, fm, 'rank', r); continue
        Y = modinv(M, p)
        def nnz(vs):
            b = kronvecs(R, vs); return int(np.count_nonzero(R.mm(b[None, :], Y)))
        print(f'{am:8s} {fm:6s} rank 256 | phys b=3/7: {nnz(phys(Fr(3,7))):3d} | generic blk1+vac: {nnz([v, u, e1, e4]):3d} | generic 2-block: {nnz([v, u, u[::-1].copy()*0+v, u]):3d}')
