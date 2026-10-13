import pickle, itertools
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
TR = list(itertools.combinations((1, 2, 3, 4), 3))
def Zp(I):
    return mp.fprod([z[i - 1] for i in I])
def w(n, I, u, j0, ths):
    """dual-Q value with the z-part: Z_I^{(u-ths)/h} q_I(u)"""
    return Zp(I) ** ((u - ths) / h) * peval(allQ[n][I], u)
# A site: tau1(thA) vs ratio of a single triple-index Q
for s in (0, 2):
    for I in TR:
        for (j1, j2) in [(-1, 0), (0, 1), (-2, -1), (-1, 1), (0, -2)]:
            Ns = [states[n].tau(1, th[s]) * w(n, I, th[s] + j2 * h, 0, th[s]) / w(n, I, th[s] + j1 * h, 0, th[s]) for n in range(256)]
            sp = max(abs(x / Ns[0] - 1) for x in Ns)
            if sp < 1e-30:
                print('A site s=%d tau1(thA) = N * QI(thA%+dh)/QI(thA%+dh), I=%s, N=%s, spread %.1e' % (s + 1, j1, j2, I, mp.nstr(Ns[0], 10), float(sp)))
# F site: tau_a(thF + a h) = N_a * det[w_l(thF + j h)]_{l in Ls, j in {0..3}\{a}} / det[...]_{j in {1,2,3}}  (Cramer of the dual Baxter)
for s in (1, 3):
    for Ls in itertools.combinations(TR, 3):
        row = []
        ok = True
        for a in (1, 2, 3):
            num_pts = [j for j in (0, 1, 2, 3) if j != a]
            Ns = []
            for n in range(256):
                Dn = mp.det(mp.matrix([[w(n, I, th[s] + j * h, 0, th[s]) for j in num_pts] for I in Ls]))
                D0 = mp.det(mp.matrix([[w(n, I, th[s] + j * h, 0, th[s]) for j in (1, 2, 3)] for I in Ls]))
                Ns.append(states[n].tau(a, th[s] + a * h) * D0 / Dn)
            sp = max(abs(x / Ns[0] - 1) for x in Ns)
            row.append('a=%d N=%s spread %.1e' % (a, mp.nstr(Ns[0], 8), float(sp)))
        print('F site s=%d triples %s: ' % (s + 1, Ls) + ' | '.join(row))
