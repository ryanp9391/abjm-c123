"""scan.py -- exact reference scan at the ABJM point z = (3/2, y, y, 4/7), y = 5/4: which product references (v (x) u) admit K_kappa?
Python fallback for reference-admissibility-Kkappa.wb sections 4, 5c.  Output: scan_out.txt"""
from rm_base import *
import random, sys
y = Fr(5, 4); zz = [Fr(3, 2), y, y, Fr(4, 7)]
zg = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]; ww = [Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)]

def tw_u(w, z):
    return [prod_([z[j] - w[m] for m in range(3)]) for j in range(4)]

def prod_(xs):
    p = Fr(1)
    for x in xs: p *= x
    return p

def tw_form(w, z):
    v = []
    for j in range(4):
        p = Fr(1)
        for m in range(4):
            if m != j: p *= (w[0] - z[m]) / (z[j] - z[m])
        v.append(p)
    u = []
    for j in range(4):
        p = Fr(1)
        for m in range(3): p *= (z[j] - w[m])
        u.append(p)
    return v, u

def comp_v(z):
    out = []
    for j in range(4):
        p = Fr(1)
        for m in range(4):
            if m != j: p *= (z[j] - z[m])
        out.append(1 / p)
    return out

def pred(v, u):
    return (v[0] * u[0] == 0 and v[3] * u[3] == 0 and v[1] * u[1] + v[2] * u[2] == 0) or (v[1] == v[2] == u[1] == u[2] == 0)

def row(name, v, u):
    v = vec(*v); u = vec(*u)
    rk, rkT, ok, cs = test(v, u, zz)
    comps = kap_comps(v, u)
    # kappa-components at class level (z2 = z3): D, B12(+B13), B14, ...: count nonzero classes
    cl = [[0], [1, 2], [3]]
    nzc = 0
    for a in range(3):
        for b in range(3):
            if any(v[i] * u[j] != 0 for i in cl[a] for j in cl[b]): nzc += 1
    dot = sum(x * y_ for x, y_ in zip(v, u))
    return (name, [str(x) for x in v], [str(x) for x in u], dot == 0, pred(v, u), ok, pred(v, u) == ok, rk, rkT)

out = open('scan_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

P('== named rows (z2 = z3 chain, pair level, exact Fractions)')
rows = [
 ('(1) W_a a=1/3', (1, 0, 0, Fr(1, 3)), (Fr(-1, 3), 0, 0, 1)),
 ('(1) W_a a=-1', (1, 0, 0, -1), (1, 0, 0, 1)),
 ('(1) bare vac', (1, 0, 0, 0), (0, 0, 0, 1)),
 ('(2) random BPS', (1, 2, -1, 3), (Fr(5, 3), Fr(-1, 6), Fr(4, 3), 0)),
 ('(2) random BPS v=(1,v2,v3,v4)', (1, Fr(-2, 3), Fr(5, 7), Fr(3, 2)), (Fr(-19, 21), Fr(1, 4), Fr(-3, 5), 1)),
 ('(3) <0|comp formal (z3=2/3)', tuple(comp_v(zg)), (1, 1, 1, 1)),
 ('(3) <0|comp z2=z3 limit', (0, 1, -1, 0), (1, 1, 1, 1)),
]
vf, uf = tw_form(ww, zg)
rows.append(('(4) <0_w| formal (z3=2/3) used on the z2=z3 chain', tuple(vf), tuple(uf)))
uL = tw_u(ww, zz)
rows.append(('(4) <0_w| z2=z3 limit generic w: (0,1,-1,0) (x) u_w', (0, 1, -1, 0), tuple(uL)))
for nm, w_ in [('w1=z1', [zz[0], Fr(9, 8), Fr(3, 7), Fr(5, 9)]), ('w1=z4', [zz[3], Fr(9, 8), Fr(3, 7), Fr(5, 9)]), ('w2=w3', [Fr(7, 5), Fr(9, 8), Fr(9, 8), Fr(5, 9)]), ('w2=w3=y', [Fr(7, 5), y, y, Fr(5, 9)]), ('w4=z4', [Fr(7, 5), Fr(9, 8), Fr(3, 7), zz[3]]), ('w1=z1,w2=w3=y', [zz[0], y, y, Fr(5, 9)])]:
    u_ = tw_u(w_, zz)
    vlim = {'w1=z1': (1, 0, 0, 0), 'w1=z4': (0, 0, 0, 1)}.get(nm, (0, 1, -1, 0))
    if nm == 'w1=z1,w2=w3=y': vlim = (1, 0, 0, 0)
    rows.append(('(4) <0_w| limit, locus ' + nm, vlim, tuple(u_)))
rows.append(('(4) <0_w| locus w1=y: v=(0,2,3,0) (alpha,beta), u2=u3=0', (0, 2, 3, 0), tuple(tw_u([y, Fr(9, 8), Fr(3, 7), Fr(5, 9)], zz))))
for r in rows:
    P(row(*r))

P('== sub-families')
sub = [
 ('Borel: W-type v=(1,0,0,2), u=(-6,0,0,3)', (1, 0, 0, 2), (-6, 0, 0, 3)),
 ('Borel: v=(1,0,0,2), u=(-6,1,2,3)', (1, 0, 0, 2), (-6, 1, 2, 3)),
 ('Borel: v generic, u=(-v4,0,0,1)', (1, Fr(-2, 3), Fr(5, 7), Fr(3, 2)), (Fr(-3, 2), 0, 0, 1)),
 ('Borel: gl2-sym v=(1,2,2,3), u=(-4,1/4,1/4,1)', (1, 2, 2, 3), (-4, Fr(1, 4), Fr(1, 4), 1)),
 ('Borel: gl2-sym v=(1,2,2,0), u=(0,0,0,1)', (1, 2, 2, 0), (0, 0, 0, 1)),
 ('Borel: v=(1,2,3,0), u=(0,3,-2,5) [p.q=0]', (1, 2, 3, 0), (0, 3, -2, 5)),
 ('Borel: v=e1, u=(0,2,3,5)', (1, 0, 0, 0), (0, 2, 3, 5)),
 ('Borel: v=(1,2,3,0), u=(-5,1,1,3) [p.q=5]', (1, 2, 3, 0), (-5, 1, 1, 3)),
 ('Borel: v=(1,2,3,2), u=(-2,3,-2,1) [p.q=0, v4u4!=0]', (1, 2, 3, 2), (-2, 3, -2, 1)),
 ('closure(u4=0): v=(1,2,3,2), u=(0,3,-2,0)', (1, 2, 3, 2), (0, 3, -2, 0)),
 ('closure(u4=0): v=(1,0,0,2), u=(0,5,7,0)', (1, 0, 0, 2), (0, 5, 7, 0)),
 ('closure(v1=0): v=(0,2,3,5), u=(7,3,-2,0)', (0, 2, 3, 5), (7, 3, -2, 0)),
 ('closure(v1=0): v=(0,2,3,0), u=(1,3,-2,4)', (0, 2, 3, 0), (1, 3, -2, 4)),
 ('closure: v=e4, u=(1,2,3,0)', (0, 0, 0, 1), (1, 2, 3, 0)),
]
for r in sub: P(row(*r))

P('== random test of the closed form')
random.seed(20261013)
def rnd(): return Fr(random.randint(1, 9) * random.choice([-1, 1]), random.randint(1, 4))
def samp_block():
    v = [Fr(0)] * 4; u = [Fr(0)] * 4
    c = random.choice([1, 2, 3])
    if c == 1: v[0] = rnd()
    if c == 2: u[0] = rnd()
    c = random.choice([1, 2, 3])
    if c == 1: v[3] = rnd()
    if c == 2: u[3] = rnd()
    p = [Fr(0), Fr(0)] if random.randint(0, 2) == 0 else [rnd(), rnd()]
    if p == [0, 0]:
        q = [Fr(0), Fr(0)] if random.randint(0, 2) == 0 else [rnd(), rnd()]
    else:
        t = Fr(0) if random.randint(0, 3) == 0 else rnd(); q = [t * p[1], -t * p[0]]
    v[1], v[2] = p; u[1], u[2] = q
    return v, u
def samp_w():
    v1, v4, t = rnd(), rnd(), rnd(); return [v1, Fr(0), Fr(0), v4], [t * v4, Fr(0), Fr(0), -t * v1]
def samp_rand():
    v = [Fr(0) if random.randint(0, 9) < 4 else rnd() for _ in range(4)]; u = [Fr(0) if random.randint(0, 9) < 4 else rnd() for _ in range(4)]
    js = [j for j in range(4) if v[j] != 0]
    if not js:
        v[random.randrange(4)] = rnd(); js = [j for j in range(4) if v[j] != 0]
    j = random.choice(js); u[j] = u[j] - sum(a * b for a, b in zip(v, u)) / v[j]
    return v, u
def tally(gen, n):
    N = pr = ac = dis = 0
    bad = []
    for _ in range(n):
        v, u = gen()
        if all(x == 0 for x in v) or all(x == 0 for x in u): continue
        assert sum(a * b for a, b in zip(v, u)) == 0
        ok = test(v, u, zz)[2]; p = pred(v, u)
        N += 1; pr += p; ac += ok
        if p != ok: dis += 1; bad.append((v, u))
    return N, pr, ac, dis, bad
for nm, g, n in [('blockwise-BPS by construction', samp_block, 60), ('W-type (p=q=0)', samp_w, 20), ('random zero pattern + v.u=0', samp_rand, 150)]:
    N, pr, ac, dis, bad = tally(g, n)
    P(nm, ': samples', N, 'predicted admissible', pr, 'K exists (exact rank test)', ac, 'disagreements', dis, bad[:3])
out.close()
