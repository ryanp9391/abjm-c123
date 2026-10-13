import numpy as np, pickle, itertools
from sympy import factorint
from core import *
from exact import *
R8 = res8(); ps = list(R8.keys())
Ys = {p: R8[p][1] for p in ps}
nz = np.argwhere(Ys[ps[0]] != 0)
th = th0; h = h0
lin_primes = set()
for i, j in itertools.combinations(range(4), 2):
    for m in range(-4, 5):
        d = th[i] - th[j] + m * h
        if d: lin_primes |= set(factorint(abs(d.numerator)).keys())
small = {2, 3, 5, 7, 11}
print('primes in theta_ij + m h (|m|<=4):', sorted(lin_primes))
vals = {}; fail = 0
for (r, c) in nz:
    v = reconstruct({p: int(Ys[p][r, c]) for p in ps})
    if v is None: fail += 1; continue
    vals[(r, c)] = v
print('entries', len(nz), 'reconstructed', len(vals), 'failed', fail)
bad = []
for k, v in vals.items():
    pr = set(factorint(abs(v.numerator)).keys()) | set(factorint(v.denominator).keys())
    extra = pr - lin_primes - small
    if extra: bad.append((k, v, extra))
print('entries with primes outside linear theta-factors:', len(bad))
for b in bad[:15]: print('  ', wordsIdx[b[0][0]], '->', wordsIdx[b[0][1]], b[1], b[2])
pickle.dump(vals, open('Y_exact.pkl', 'wb'))
