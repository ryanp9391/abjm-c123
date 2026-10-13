import itertools, pickle
import numpy as np
from numpy.polynomial import chebyshev as C
exec(open('bop.py').read().split("W = pickle.load")[0])
W = pickle.load(open('words_exact.pkl', 'rb'))['words']
Wf = np.array([[float(x) for x in row] for row in W])
deg = 24
a, b = -3.0, 3.0
xs = np.cos(np.pi * (np.arange(deg + 1) + 0.5) / (deg + 1))
us = (a + b) / 2 + (b - a) / 2 * xs
vals = np.zeros((256, deg + 1))
res = 0
for m, u in enumerate(us):
    Bu = Bop(u)
    bB = Wf @ Bu
    for n in range(256):
        p = np.argmax(np.abs(Wf[n]))
        lam = bB[n, p] / Wf[n, p]
        res = max(res, np.abs(bB[n] - lam * Wf[n]).max() / np.abs(bB[n]).max())
        vals[n, m] = lam
print('max B-eigen residual', res)
roots = {}
for n, w in enumerate(WORDS):
    c = C.chebfit(xs, vals[n], deg)
    r = C.chebroots(c)
    roots[w] = (a + b) / 2 + (b - a) / 2 * r
def ident(r):
    out = []
    for x in r:
        best = min(((abs(x - (thf[s] + j * hf)), s + 1, j) for s in range(4) for j in range(-8, 9)))
        out.append((best[1], best[2], best[0]))
    return out
for w in WORDS:
    if (w[1] == 0 and w[3] == 0 and w[2] == 0) or (w[0] == 0 and w[2] == 0 and w[3] == 0) or w in [(1,1,0,0),(3,2,1,1)]:
        idr = ident(roots[w])
        per = {s: sorted(j for (ss, j, e) in idr if ss == s) for s in range(1, 5)}
        print(w, per, 'maxerr %.1e' % max(e for (_, _, e) in idr))
pickle.dump(roots, open('Broots.pkl', 'wb'))
