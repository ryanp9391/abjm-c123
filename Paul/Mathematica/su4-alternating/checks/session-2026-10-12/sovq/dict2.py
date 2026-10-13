import pickle, itertools
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
def Qt(n, m, thA, j):
    """regularised Q_m(thA - j h), j = 1..4 (common factor z_m^{thA/h} Phi4(thA-4h) dropped)"""
    return z[m - 1] ** (-j) * mp.fprod([V(4, thA - b * h) for b in range(j, 4)]) * peval(allQ[n][(m,)], thA - j * h)
def D(n, I, thA, a):
    pts = [j for j in (1, 2, 3, 4) if j != a + 1]
    return mp.det(mp.matrix([[Qt(n, m, thA, j) for j in pts] for m in I]))
for site, thA in [('A1', th[0]), ('A3', th[2])]:
    for I in itertools.combinations((1, 2, 3, 4), 3):
        row = []
        for a in (1, 2, 3):
            Ns = [states[n].tau(a, thA - h) * D(n, I, thA, 0) / D(n, I, thA, a) for n in range(256)]
            spread = max(abs(x - Ns[0]) for x in Ns) / abs(Ns[0])
            row.append('a=%d N=%s spread %.1e' % (a, mp.nstr(Ns[0], 12), float(spread)))
        print(site, I, ' | '.join(row))
