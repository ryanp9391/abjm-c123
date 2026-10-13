import itertools, pickle
import numpy as np
from sov_l4 import *
hf = 0.5; thf = tuple(float(x) for x in TH0); zf = tuple(float(x) for x in Z0)
def Tm(u):
    P = monodromy_dicts(u, thf, hf)
    return [[to_dense(P[i][j], float) for j in range(4)] for i in range(4)]
Tc = {}
def T(u):
    if u not in Tc: Tc[u] = Tm(u)
    return Tc[u]
def sgn(p):
    s = 1
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s
def qmin(rows, cols, v):
    a = len(cols)
    tot = np.zeros((256, 256))
    for p in itertools.permutations(range(a)):
        M = np.identity(256)
        for k in range(a):
            M = M @ T(v - k * hf)[rows[p[k]]][cols[k]]
        tot += sgn(p) * M
    return tot
def Bop(u):
    tot = np.zeros((256, 256))
    for k1 in (1, 2, 3):
        for K2 in itertools.combinations((1, 2, 3), 2):
            A = T(u)[k1 - 1][0]
            Bm = qmin([k - 1 for k in K2], [0, k1], u + hf)
            C = qmin([0, 1, 2], [0, K2[0], K2[1]], u + 2 * hf)
            tot += A @ Bm @ C
    return tot
W = pickle.load(open('words_exact.pkl', 'rb'))['words']
Wf = np.array([[float(x) for x in row] for row in W])
us = [0.1 * k + 0.0371 for k in range(-6, 7)]  # 13 points for degree 12
vals = np.zeros((256, len(us)))
res = 0
for m, u in enumerate(us):
    Bu = Bop(u)
    for n in range(256):
        b = Wf[n]; bB = b @ Bu
        p = np.argmax(np.abs(b))
        lam = bB[p] / b[p]
        res = max(res, np.abs(bB - lam * b).max() / np.abs(bB).max())
        vals[n, m] = lam
print('max B-eigen residual of words', res)
Vm = np.array([[u ** k for k in range(13)] for u in us])
roots = {}
for n, w in enumerate(WORDS):
    c = np.linalg.solve(Vm, vals[n])
    r = np.roots(c[::-1])
    roots[w] = np.sort_complex(r)
# identify roots as theta_s + j h
def ident(r):
    out = []
    for x in r:
        best = min(((abs(x - (thf[s] + j * hf)), s + 1, j) for s in range(4) for j in range(-6, 7)))
        out.append((best[1], best[2], best[0]))
    return out
for w in WORDS:
    idr = ident(roots[w])
    if w[1] == 0 and w[3] == 0 or (w[0] == 0 and w[2] == 0):
        per = {s: sorted(j for (ss, j, e) in idr if ss == s) for s in range(1, 5)}
        print(w, per, 'maxerr %.1e' % max(e for (_, _, e) in idr))
pickle.dump(roots, open('Broots.pkl', 'wb'))
