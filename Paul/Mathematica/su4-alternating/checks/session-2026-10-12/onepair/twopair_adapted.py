import numpy as np, pickle, random
from fractions import Fraction
from itertools import combinations
from fastchain import *
d = pickle.load(open('twopair_W.pkl', 'rb')); labs = d['labs']
Winv = pickle.load(open('twopair_Winv.pkl', 'rb'))
zA = [Fraction(3, 2), Fraction(5, 4), Fraction(2, 3), Fraction(4, 7)]
Rm = Ring('mod')
c = [Rm.c(sum(np.prod(s) for s in combinations(zA, k))) for k in range(1, 5)]
G = np.array([[c[0], (-c[1]) % P, c[2], (-c[3]) % P], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]], dtype=object)
GT = G.T.copy()
def matpoly(M, coeffs):
    out = np.zeros((4, 4), dtype=object); Mk = np.eye(4, dtype=object)
    for a in coeffs: out = (out + a * Mk) % P; Mk = Mk.dot(M) % P
    return out
e1 = np.array([1, 0, 0, 0], dtype=object); e4 = np.array([0, 0, 0, 1], dtype=object)
def coeffs(v1, u1, v2, u2):
    x = np.kron(np.kron(np.kron(v1, u1), v2), u2) % P
    cc = x.dot(Winv) % P
    C = np.zeros((16, 16), dtype=object)
    for n, (a1, k2, a3, k4) in enumerate(labs): C[4 * a1 + k2, 4 * a3 + k4] = cc[n]
    return C
random.seed(11)
r = lambda: random.randrange(1, P)
def A(deg): return e1.dot(matpoly(GT, [r() for _ in range(deg + 1)])) % P
def F(deg): return e4.dot(matpoly(G, [r() for _ in range(deg + 1)])) % P
cases = [('p,q deg (1,0) | vac', (A(1), F(0), A(0), F(0))),
         ('p,q deg (0,1) | vac', (A(0), F(1), A(0), F(0))),
         ('p,q deg (1,1) | vac', (A(1), F(1), A(0), F(0))),
         ('vac | p,q deg (1,0)', (A(0), F(0), A(1), F(0))),
         ('vac | p,q deg (0,1)', (A(0), F(0), A(0), F(1))),
         ('vac | p,q deg (1,1)', (A(0), F(0), A(1), F(1))),
         ('(1,1) | (1,1)', (A(1), F(1), A(1), F(1))),
         ('(1,0) | (1,0)', (A(1), F(0), A(1), F(0))),
         ('(0,1) | (0,1)', (A(0), F(1), A(0), F(1))),
         ('(2,2) | (2,2)', (A(2), F(2), A(2), F(2)))]
for name, vs in cases:
    C = coeffs(*vs)
    sup = [labs[n] for n in range(256) if C[4 * labs[n][0] + labs[n][1], 4 * labs[n][2] + labs[n][3]] % P]
    print('%-22s rank %2d, #words %3d' % (name, rank_mod(C), len(sup)), sup if len(sup) <= 14 else '')
