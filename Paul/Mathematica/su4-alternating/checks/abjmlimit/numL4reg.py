"""delta = 0 (pair-by-pair ABJM point): regularised basis overlaps on eigenstates; residue comparison."""
import pickle, time
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
from words import WIDX
from numL4 import to_mp, to_np, refine, B_TEST

mp.mp.dps = 60
t0 = time.time()
R = pickle.load(open('exact_R0.pkl', 'rb'))
T = pickle.load(open('exact_t1_reg.pkl', 'rb')) if False else pickle.load(open('exact_t1_0.pkl', 'rb'))
Rmp = to_mp(R)
bmp = mp.matrix([mp.mpf(x.numerator) / x.denominator for x in B_TEST])
r = mp.lu_solve(Rmp.T, bmp)
Tmp = to_mp(T)
Tnp = to_np(T)
ev, evec = np.linalg.eig(Tnp)
gaps = np.array([np.min(np.abs(np.delete(ev, i) - ev[i])) for i in range(len(ev))])
order = np.argsort(-gaps)[:8]
I = lambda w: WIDX.index(w)
rows = []
for i in order:
    lam, psi, res = refine(Tmp, Tnp, ev[i], evec[:, i])
    ov = Rmp * psi
    o0 = ov[I((0, 0, 0, 0))]
    tF2 = ov[I((0, 1, 0, 0))] / o0
    tF4 = ov[I((0, 0, 0, 1))] / o0
    x1, w1, y1 = [ov[I((a, 0, 0, 0))] / o0 for a in (1, 2, 3)]
    x3, w3, y3 = [ov[I((0, 0, a, 0))] / o0 for a in (1, 2, 3)]

    def lamA(x, w, y, tF, a, k):
        return [1, x, (w if k <= 2 else x * tF), y][a]

    pred = [lamA(x1, w1, y1, tF2, a1, k2) * tF2 ** k2 * lamA(x3, w3, y3, tF4, a3, k4) * tF4 ** k4 for (a1, k2, a3, k4) in WIDX]
    fact = max(abs(ov[j] / o0 - pred[j]) / max(1, abs(pred[j])) for j in range(256))
    S = mp.fsum(r[j] * pred[j] for j in range(256))
    direct = mp.fsum(bmp[j] * psi[j] for j in range(256)) / o0
    rows.append((mp.nstr(res, 3), mp.nstr(fact, 3), mp.nstr(S, 12), mp.nstr(direct, 12), mp.nstr(abs(S - direct) / abs(direct), 3),
                 mp.nstr(x1, 8), mp.nstr(w1, 8), mp.nstr(y1, 8), mp.nstr(tF2, 8)))
print('d = 2h (pair-by-pair), regularised basis: max|r| =', mp.nstr(max(abs(v) for v in r), 6), ' time', round(time.time() - t0, 1))
for rr in rows:
    print('  eig-res %s factorisation %s  sum r*lambda = %s  direct = %s  rel.diff %s | x1=%s w1=%s y1=%s tauF2=%s' % rr)
# residue comparison with delta runs
out = {'r': [str(v) for v in r], 'rows': rows}
for tag in ('100', '1000', '10000'):
    try:
        d = pickle.load(open('numL4_%s.pkl' % tag, 'rb'))
    except Exception:
        continue
    c = [mp.mpf(s) if 'j' not in s else mp.mpc(complex(s.replace('(', '').replace(')', ''))) for s in d['c']]
    delta = mp.mpf(1) / int(tag)
    sel = [(2, 0, 2, 0), (3, 0, 3, 0), (3, 1, 2, 2), (2, 1, 0, 0), (0, 2, 3, 3), (0, 0, 0, 0), (1, 0, 2, 3)]
    vals = []
    for w in sel:
        n = (w[0] == 3 or (w[0] == 2 and w[1] <= 2)) + (w[2] == 3 or (w[2] == 2 and w[3] <= 2))
        vals.append((w, n, mp.nstr(delta ** n * c[I(w)], 10), mp.nstr(r[I(w)], 10)))
    # finite combination: c(2,k,0,0)+c(1,k+1,0,0) -> r(1,k+1,0,0)
    comb = [(k, mp.nstr(c[I((2, k, 0, 0))] + c[I((1, k + 1, 0, 0))], 10), mp.nstr(r[I((1, k + 1, 0, 0))], 10)) for k in range(3)]
    print('delta = 1/%s: max|c|*delta^2 = %s' % (tag, mp.nstr(max(abs(v) for v in c) * delta ** 2, 6)))
    for v in vals:
        print('   word', v[0], 'n =', v[1], ' delta^n c =', v[2], '  r =', v[3])
    for v in comb:
        print('   c(2,%d,0,0)+c(1,%d,0,0) = %s   r(1,%d,0,0) = %s' % (v[0], v[0] + 1, v[1], v[0] + 1, v[2]))
pickle.dump(out, open('numL4reg.pkl', 'wb'))
