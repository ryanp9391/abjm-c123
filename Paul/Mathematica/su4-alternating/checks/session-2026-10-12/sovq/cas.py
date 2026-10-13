import pickle, itertools
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
qs = pickle.load(open('qs.pkl', 'rb'))
def Pcas(n, I):
    """P_I(u) = det_{s,m}[ z_m^{|I|-1-s} q_m(u - s h) ], m in I (increasing), s = 0..|I|-1"""
    k = len(I)
    M = [[pscale(pshift(qs[n][m][0], -s * h), z[m - 1] ** (k - 1 - s)) for m in I] for s in range(k)]
    return ptrim(pdet(M))
vac = 255
for k in (2, 3, 4):
    for I in itertools.combinations(range(1, 5), k):
        P = Pcas(vac, I)
        rs = proots(P)
        print(I, 'deg', len(P) - 1, 'roots', [mp.nstr(r, 8) for r in rs])
print('theta', [mp.nstr(t, 8) for t in th])
