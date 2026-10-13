import pickle, itertools
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
TR = list(itertools.combinations((1, 2, 3, 4), 3))
def Zp(I): return mp.fprod([z[i - 1] for i in I])
def w(n, I, u, ths): return Zp(I) ** ((u - ths) / h) * peval(allQ[n][I], u)
Ls = TR[:3]
sample = range(0, 256, 9)
for s in (1, 3):
    for c in range(-4, 3):
        row = []
        for a in (1, 2, 3):
            Ns = []
            for n in sample:
                M = [[w(n, I, th[s] + (c + j) * h, th[s]) for j in range(4)] for I in Ls]
                Dn = mp.det(mp.matrix([[M[l][j] for j in range(4) if j != a] for l in range(3)]))
                D0 = mp.det(mp.matrix([[M[l][j] for j in (1, 2, 3)] for l in range(3)]))
                Ns.append(states[n].tau(a, th[s] + a * h) * D0 / Dn)
            sp = max(abs(x / Ns[0] - 1) for x in Ns)
            row.append('a=%d spread %.1e' % (a, float(sp)))
        print('F s=%d offset c=%d: ' % (s + 1, c) + ' | '.join(row))
