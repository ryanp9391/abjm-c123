import os, tempfile
import numpy as np, io, contextlib, collections, pickle, json, itertools
with contextlib.redirect_stdout(io.StringIO()):
    from inv import *
    from symm import coeffs, res
    from polyverify import Fpoly, reduce_poly
from core import S

THR = 1e-9
NAMES = ['a2', 'a3', 'd2', 'd3', 'mv', 'mu', 'p22', 'p23', 'p32', 'p33']
OUT = {}


def num(x):
    if abs(x - round(x)) < 1e-9:
        return str(int(round(x)))
    return '%.6g' % x


def cstr(c):
    r = c.real if abs(c.real) > 1e-11 else 0.0
    i = c.imag if abs(c.imag) > 1e-11 else 0.0
    if i == 0:
        return num(r)
    if r == 0:
        return num(i) + ' I'
    return '(' + num(r) + (' - ' if i < 0 else ' + ') + num(abs(i)) + ' I)'


def monstr(key):
    parts = []
    for n, e in zip(NAMES, key):
        if e == 1:
            parts.append(n)
        elif e > 1:
            parts.append('%s^%d' % (n, e))
    return ' '.join(parts)


def redmon(key):
    a2, a3, d2, d3, mv, mu, p22, p23, p32, p33 = key
    parts = []
    for n, e in (('x', p23), ('y', p32), ('μv', mv), ('μu', mu)):
        if e == 1:
            parts.append(n)
        elif e > 1:
            parts.append('%s^%d' % (n, e))
    return ' '.join(parts)


def nz(L, k, M, bps):
    c = coeffs(L, k, M, bps)
    return [(key, v) for key, v in c.items() if abs(v) > THR]


def formstr(L, k, M, bps, red=False):
    t = nz(L, k, M, bps)
    if not t:
        return '0'
    parts = []
    for key, v in t:
        m = redmon(key) if red else monstr(key)
        parts.append(cstr(v) + (' ' + m if m else ''))
    return '  +  '.join(parts)


def dynstr(L, k):
    return ','.join(str(int(x)) for x in S[L][k]['dyn'])


def Estr(L, k):
    return '%.4f' % S[L][k]['E']


def table_md(L, Ms, bpsboth=False):
    rows = []
    if bpsboth:
        hdr = '| # | Dynkin | E_MZ | M | F, BPS pairs (p22 = p33 = 0) | F, no BPS condition |\n|---|---|---|---|---|---|'
    else:
        hdr = '| # | Dynkin [p,q,r] | E_MZ | M | F (BPS pairs) in the invariants | F/F_vac in x, y, μv, μu |\n|---|---|---|---|---|---|'
    for k in range(len(S[L])):
        for M in Ms:
            if bpsboth:
                if nz(L, k, M, False):
                    rows.append('| %d | %s | %s | %d | %s | %s |' % (k + 1, dynstr(L, k), Estr(L, k), M, formstr(L, k, M, True), formstr(L, k, M, False)))
            else:
                if nz(L, k, M, True):
                    rows.append('| %d | %s | %s | %d | %s | %s |' % (k + 1, dynstr(L, k), Estr(L, k), M, formstr(L, k, M, True), formstr(L, k, M, True, True)))
    return hdr + '\n' + '\n'.join(rows), len(rows)


OUT['table2'], n2 = table_md(2, range(3), True)
OUT['table3'], n3 = table_md(3, range(4))
OUT['table4'], n4 = table_md(4, range(5))
OUT['nrows'] = (n2, n3, n4)

# polynomials: counts
cnt = {}
for L in (2, 3):
    mx = 0
    zero = 0
    tot = 0
    for k in range(len(S[L])):
        for M in range(L + 1):
            F = Fpoly(S[L][k]['psi'], L, M)
            nt = len([1 for c in F.values() if abs(c) > 1e-12])
            mx = max(mx, nt)
            zero += (nt == 0)
            tot += 1
    cnt[L] = (tot, mx, zero)
OUT['polys'] = cnt

# fit summary
fs = {}
for L in (2, 3, 4):
    ent = [(key, r) for key, r in res.items() if key[0] == L]
    withm = [r for key, r in ent if r['nm']]
    without = [r for key, r in ent if not r['nm']]
    fs[L] = (len(ent), len(withm), max(r['res'] for r in withm), len(without), max(r['Fmax'] for r in without),
             sum(1 for r in withm if r['rank'] < r['nm']))
OUT['fits'] = fs


# structure counts and classes
def nstruct(L, M, p, q, r):
    n2, n3 = L - M, M
    if (2 * q + r - p) % 4 or (p + 2 * q - r) % 4:
        return 0
    e5 = (2 * q + r - p) // 4
    e6 = (p + 2 * q - r) // 4
    s2 = p + 2 * q + r
    if e5 < 0 or e6 < 0 or s2 % 2:
        return 0
    s = L - s2 // 2
    if s < 0:
        return 0
    return max(0, min(s, n2 - e5, n3 - e6) - max(0, s - n3 + e5, s - n2 + e6) + 1)


rows = []
bad = 0
for L in (2, 3, 4):
    dyns = sorted(set(tuple(int(x) for x in s['dyn']) for s in S[L]))
    for d in dyns:
        p, q, r = d
        cntd = sum(1 for s in S[L] if tuple(int(x) for x in s['dyn']) == d)
        e5 = (2 * q + r - p) / 4
        e6 = (p + 2 * q - r) / 4
        s_ = L - (p + 2 * q + r) / 2
        Ns = [nstruct(L, M, p, q, r) for M in range(L + 1)]
        rows.append('| %d | %d,%d,%d | %d | %s | %s | %s | %s |' % (L, p, q, r, cntd, num(e5), num(e6), num(s_), ' '.join(map(str, Ns))))
for L in (2, 3, 4):
    for k, s in enumerate(S[L]):
        p, q, r = [int(x) for x in s['dyn']]
        for M in range(L + 1):
            if nstruct(L, M, p, q, r) != res[(L, k, M, True)]['nm']:
                bad += 1
OUT['dyntab'] = '| L | [p,q,r] | # states | e5 (power of mv) | e6 (power of mu) | s = e8+e9 (total pairing power) | # invariant structures N for M = 0..L |\n|---|---|---|---|---|---|---|\n' + '\n'.join(rows)
OUT['nstruct_bad'] = bad


def classify(L, k, M):
    r = res[(L, k, M, True)]
    if r['nm'] == 0:
        return 'kinematic zero'
    t = nz(L, k, M, True)
    if not t:
        return 'dynamical zero'
    supp = [(key[7], key[8]) for key, v in t]
    if sum(supp[0]) == 0:
        return 'no pairing (s=0)'
    e8 = min(a for a, b in supp)
    e9 = min(b for a, b in supp)
    if e8 >= 1 and e9 >= 1:
        return 'needs both p23,p32'
    if e8 == 0 and e9 >= 1:
        return 'survives p23=0 (needs p32)'
    if e9 == 0 and e8 >= 1:
        return 'survives p32=0 (needs p23)'
    return 'survives either pairing = 0'


CL = ['kinematic zero', 'dynamical zero', 'no pairing (s=0)', 'needs both p23,p32', 'survives p23=0 (needs p32)', 'survives p32=0 (needs p23)', 'survives either pairing = 0']
crow = []
dz = []
for L in (2, 3, 4):
    for M in range(L + 1):
        c = collections.Counter(classify(L, k, M) for k in range(len(S[L])))
        crow.append('| %d | %d | ' % (L, M) + ' | '.join(str(c.get(x, 0)) for x in CL) + ' |')
    for k, s in enumerate(S[L]):
        for M in range(L + 1):
            if classify(L, k, M) == 'dynamical zero':
                dz.append((L, k + 1, dynstr(L, k), Estr(L, k), M, res[(L, k, M, True)]['nm']))
OUT['classtab'] = '| L | M | ' + ' | '.join(CL) + ' |\n|' + '---|' * (len(CL) + 2) + '\n' + '\n'.join(crow)
OUT['dz'] = dz

# checks: note polarisation, cyclicity, x<->y symmetry
e = np.eye(4)
note = [e[0], e[3], e[0] - e[3], e[3] + e[0]]
worst = 0
for L in (2, 3, 4):
    for k, s in enumerate(S[L]):
        for M in range(L + 1):
            worst = max(worst, abs(Fval(s['psi'], L, M, note) - s['C'][2 * (L - M)]))
OUT['note_diff'] = worst
rng2 = np.random.default_rng(5)
worstc = 0
for L in (2, 3, 4):
    for k, s in enumerate(S[L]):
        for M in range(L + 1):
            pt = [x / np.linalg.norm(x) for x in [rng2.normal(size=4) + 1j * rng2.normal(size=4) for _ in range(4)]]
            a = Fval(s['psi'], L, M, pt)
            b = Fval(s['psi'], L, L - M, [pt[2], pt[3], pt[0], pt[1]])
            worstc = max(worstc, abs(a - b))
OUT['cyc_diff'] = worstc


def swapQ(m):
    a2, a3, d2, d3, mv, mu, p22, p23, p32, p33 = m
    return (d2, d3, a2, a3, mu, mv, p22, p32, p23, p33)


worstq = 0
ncase = 0
for L in (2, 3, 4):
    for k, s in enumerate(S[L]):
        p, q, r = [int(x) for x in s['dyn']]
        if p != r:
            continue
        for M in range(L + 1):
            c = coeffs(L, k, M, True)
            for key, v in c.items():
                worstq = max(worstq, abs(v - c.get(swapQ(key), 0)))
            if c:
                ncase += 1
OUT['sym_diff'] = (worstq, ncase)

# visible polarisations
P = [np.array(x, dtype=complex) for x in ([1, 0, 0, 0], [0, 0, 0, 1], [1, 1, 0, 1], [1, 1, 0, -2])]
Pn = note
vis = collections.defaultdict(list)
worstF = 0
for L in (2, 3, 4):
    for k, s in enumerate(S[L]):
        for M in range(1, L):
            Iv = invs(*P)
            Fv = (Iv['a2'] * Iv['d2']) ** (L - M) * (Iv['a3'] * Iv['d3']) ** M
            Fs = Fval(s['psi'], L, M, P) / Fv
            In = invs(*Pn)
            Fvn = (In['a2'] * In['d2']) ** (L - M) * (In['a3'] * In['d3']) ** M
            Fn = Fval(s['psi'], L, M, Pn) / Fvn
            # invariant form at P*
            tot = 0
            for key, v in nz(L, k, M, True):
                tot += v * np.prod([Iv[n] ** e_ for n, e_ in zip(NAMES, key)]) / Fv
            worstF = max(worstF, abs(tot - Fs))
            if abs(Fn) < 1e-9 and abs(Fs) > 1e-9:
                vis[(L, dynstr(L, k), M)].append(Fs)
OUT['Fhat_check'] = worstF
vrows = []
for (L, d, M), lst in sorted(vis.items()):
    vrows.append('| %d | %s | %d | %d | %s |' % (L, d, M, len(lst), cstr(lst[0])))
OUT['vis'] = '| L | [p,q,r] | M | # states with F/F_vac(W-1) = 0 but F/F_vac(P*) != 0 | example F/F_vac(P*) |\n|---|---|---|---|---|\n' + '\n'.join(vrows)

# x,y reduction: classes statistics
st = collections.Counter()
for L in (2, 3, 4):
    for k, s in enumerate(S[L]):
        for M in range(L + 1):
            t = nz(L, k, M, True)
            if t:
                st[(L, len(t))] += 1
OUT['nstruct_nonzero'] = dict(st)
json.dump({k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()}) for k, v in OUT.items()}, open(os.path.join(tempfile.gettempdir(), 'twoblock_nbdata.json'), 'w'), default=str)
for k in ('nrows', 'polys', 'fits', 'nstruct_bad', 'note_diff', 'cyc_diff', 'sym_diff', 'Fhat_check', 'nstruct_nonzero'):
    print(k, OUT[k])
print(OUT['classtab'])
print(OUT['vis'])
print(len(OUT['dz']))
