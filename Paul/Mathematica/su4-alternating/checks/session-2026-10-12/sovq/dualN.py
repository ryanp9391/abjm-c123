"""State-independent constants of the ket (dual) SoV wave functions, dual Q's taken as w_I(u) = Z_I^{(u-th_s)/h} q_I(u):
  A site s:  tau1(th_s)        = NA_s   * w_I(th_s - 2h)/w_I(th_s - h)                          (any triple I)
  F site s:  tau_n(th_s + n h) = NF_s,n * det[w_I(th_s+(j-1)h)]_{j in {0..3}\{n}} / det[...]_{j in {1,2,3}}  (any 3 triples)"""
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
NA = {}; NF = {}
for s in (0, 2):
    n = 5
    NA[s] = states[n].tau(1, th[s]) * w(n, TR[0], th[s] - h, th[s]) / w(n, TR[0], th[s] - 2 * h, th[s])
for s in (1, 3):
    for a in (1, 2, 3):
        n = 5
        Ls = TR[:3]
        Mx = [[w(n, I, th[s] + (j - 1) * h, th[s]) for j in range(4)] for I in Ls]
        Dn = mp.det(mp.matrix([[Mx[l][j] for j in range(4) if j != a] for l in range(3)]))
        D0 = mp.det(mp.matrix([[Mx[l][j] for j in (1, 2, 3)] for l in range(3)]))
        NF[(s, a)] = states[n].tau(a, th[s] + a * h) * D0 / Dn
for k, v in NA.items(): print('NA site', k + 1, mp.nstr(v, 20))
for k, v in NF.items(): print('NF site', k[0] + 1, 'n =', k[1], mp.nstr(v, 20), mp.identify(v))
pickle.dump({'NA': NA, 'NF': NF}, open('dualN.pkl', 'wb'))
