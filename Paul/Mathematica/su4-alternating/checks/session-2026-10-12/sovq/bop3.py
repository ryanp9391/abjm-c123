import itertools, pickle
import numpy as np
exec(open('bop.py').read().split("W = pickle.load")[0])
W = pickle.load(open('words_exact.pkl', 'rb'))['words']
Wf = np.array([[float(x) for x in row] for row in W])
us = np.linspace(-2.9, 2.9, 41) + 0.0123
vals = np.zeros((256, len(us)))
for m, u in enumerate(us):
    bB = Wf @ Bop(u)
    for n in range(256):
        p = np.argmax(np.abs(Wf[n]))
        vals[n, m] = bB[n, p] / Wf[n, p]
pickle.dump((us, vals), open('Bvals.pkl', 'wb'))
def ratfit(n, d):
    r = vals[n] / vals[0]
    # r(u) * Q(u) - P(u) = 0, P,Q degree d, Q monic
    A = np.array([[r[m] * us[m] ** k for k in range(d + 1)] + [-us[m] ** k for k in range(d + 1)] for m in range(len(us))])
    U, S, Vt = np.linalg.svd(A)
    v = Vt[-1]
    Q = v[:d + 1]; P = v[d + 1:]
    return S[-1] / S[0], np.roots(P[::-1]), np.roots(Q[::-1]), P[-1] / Q[-1]
def ident(x):
    best = min(((abs(x - (thf[s] + j * hf)), s + 1, j) for s in range(4) for j in range(-8, 9)))
    return 'th%d%+dh(%.0e)' % (best[1], best[2], best[0])
for w in [(1,0,0,0),(2,0,0,0),(3,0,0,0),(0,1,0,0),(0,2,0,0),(0,3,0,0),(0,0,1,0),(0,0,0,1),(1,1,0,0),(3,3,3,3)]:
    n = WORDS.index(w)
    for d in range(1, 8):
        s, zr, pr, c = ratfit(n, d)
        if s < 1e-9:
            break
    # cancel common roots
    zr = list(zr); pr = list(pr)
    for x in list(zr):
        for y in list(pr):
            if abs(x - y) < 1e-6 and x in zr and y in pr:
                zr.remove(x); pr.remove(y); break
    print(w, 'd=%d sv=%.1e' % (d, s), 'zeros:', [ident(x) for x in zr], 'poles:', [ident(x) for x in pr])
