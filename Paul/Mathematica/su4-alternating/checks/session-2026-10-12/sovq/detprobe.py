"""Bilinear SoV form B(F,G) = sum_{x,y} G(y) Mt_{y,x} F(x) with the Q-normalised measure, and a GLMR-type determinant probe."""
import pickle, itertools, random
import mpmath as mp
from qfun import *
from poly import *
from sov_l4 import WORDS
mp.mp.dps = 40
ms = pickle.load(open('measure_exact.pkl', 'rb')); M = ms['M']
dn = pickle.load(open('dualN.pkl', 'rb')); NA, NF = dn['NA'], dn['NF']
tomp = lambda x: mp.mpf(x.numerator) / x.denominator
def Ny(y):
    return NA[0] ** y[0] * (NF[(1, y[1])] if y[1] else 1) * NA[2] ** y[2] * (NF[(3, y[3])] if y[3] else 1)
Mt = {}
for (xp, x), v in M.items():
    y = tuple(3 - c for c in xp)
    Mt[(y, x)] = Ny(y) * tomp(v)
pickle.dump(Mt, open('Mt.pkl', 'wb'))
# site functions. bra: F site (s=1,3): fF[s] = (f(th-h), f(th)) -> F_F(k) = f0^(3-k) f1^k
#                  A site (s=0,2): fA[s][i][j] = f_i(th - j h), i=0..2, j=1..4 -> F_A(a) = det_{i, j in {1..4}\{a+1}}
#  ket: A site: gA[s] = (g(th-h), g(th-2h)) -> G_A(m) = g0^(3-m) g1^m
#       F site: gF[s][i][j] = g_i(th + (j-1)h), j = 0..3 -> G_F(n) = det_{i, j in {0..3}\{n}}
def Fx(fF, fA, x):
    p = mp.mpf(1)
    for s in (0, 2):
        cols = [j for j in (1, 2, 3, 4) if j != x[s] + 1]
        p *= mp.det(mp.matrix([[fA[s][i][j] for j in cols] for i in range(3)]))
    for s in (1, 3):
        p *= fF[s][0] ** (3 - x[s]) * fF[s][1] ** x[s]
    return p
def Gy(gA, gF, y):
    p = mp.mpf(1)
    for s in (0, 2):
        p *= gA[s][0] ** (3 - y[s]) * gA[s][1] ** y[s]
    for s in (1, 3):
        cols = [j for j in (0, 1, 2, 3) if j != y[s]]
        p *= mp.det(mp.matrix([[gF[s][i][j] for j in cols] for i in range(3)]))
    return p
def Bform(fF, fA, gA, gF):
    Fv = {x: Fx(fF, fA, x) for x in WORDS}
    Gv = {y: Gy(gA, gF, y) for y in WORDS}
    return mp.fsum(Gv[y] * v * Fv[x] for (y, x), v in Mt.items())
if __name__ == '__main__':
    # eigenstate check: B(Psi) = <PsiL|Psi> * Psi(x0) PsiT(y0) / (<0|Psi><PsiL|0bar>)
    data = pickle.load(open('eigen.pkl', 'rb')); states = [State(d) for d in data['states']]
    allQ = pickle.load(open('allQ.pkl', 'rb'))
    TR = list(itertools.combinations((1, 2, 3, 4), 3))
    J = 5
    def Qt(n, i, ths, j): return z[i - 1] ** (-j) * mp.fprod([V(4, ths - b * h) for b in range(j, J)]) * peval(allQ[n][(i,)], ths - j * h)
    def w(n, I, u, ths): return mp.fprod([z[i - 1] for i in I]) ** ((u - ths) / h) * peval(allQ[n][I], u)
    def eigdata(n, trip=(1, 2, 3), jF=1, IA=TR[0], LF=TR[:3]):
        fA = {s: [[None] + [Qt(n, i, th[s], j) for j in (1, 2, 3, 4)] for i in trip] for s in (0, 2)}
        fF = {s: (Qt(n, jF, th[s], 1), Qt(n, jF, th[s], 0)) for s in (1, 3)}
        gA = {s: (w(n, IA, th[s] - h, th[s]), w(n, IA, th[s] - 2 * h, th[s])) for s in (0, 2)}
        gF = {s: [[w(n, I, th[s] + (j - 1) * h, th[s]) for j in range(4)] for I in LF] for s in (1, 3)}
        return fF, fA, gA, gF
    lr = pickle.load(open('lr_vecs.pkl', 'rb'))
    out = []
    for n in [0, 17, 100, 255]:
        fF, fA, gA, gF = eigdata(n)
        # direct: tau-products
        Bv = Bform(fF, fA, gA, gF) / (Fx(fF, fA, (0, 0, 0, 0)) * Gy(gA, gF, (0, 0, 0, 0)))
        out.append((n, Bv))
    pickle.dump(out, open('Beig.pkl', 'wb'))
    for n, b in out: print('state', n, 'B/(Psi(x0)PsiT(y0)) =', mp.nstr(b, 15))
