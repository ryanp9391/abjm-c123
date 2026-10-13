import sympy as sp, time
from itertools import permutations, combinations
from chain import perm_sign, companion
from fractions import Fraction as Fr
d, u, eta = sp.symbols('d u eta')
h = sp.Rational(1,2)
Z = [sp.Rational(3,2), sp.Rational(5,4), sp.Rational(2,3), sp.Rational(4,7)]
Gc = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in companion([Fr(3,2),Fr(5,4),Fr(2,3),Fr(4,7)])])
def E(i,j):
    m = sp.zeros(4,4); m[i,j]=1; return m
def rho(ty,i,j): return E(i,j) if ty=='A' else -E(j,i)
types=['A','F']
def T(uu, th):
    # T_ij = sum_k L1_ik (x) L2_kj ; then twist
    L = [[[ (uu-th[s])*(1 if i==j else 0)*sp.eye(4) + h*rho(types[s],i,j) for j in range(4)] for i in range(4)] for s in range(2)]
    Tm = [[sum((sp.kronecker_product(L[0][i][k], L[1][k][j]) for k in range(4)), sp.zeros(16,16)) for j in range(4)] for i in range(4)]
    TG = [[sum((Tm[i][k]*Gc[k,j] for k in range(4)), sp.zeros(16,16)) for j in range(4)] for i in range(4)]
    return TG
def t(a, uu, th):
    if a==0: return sp.eye(16)
    Ts = [T(uu - k*h, th) for k in range(a)]
    acc = sp.zeros(16,16)
    for I in combinations(range(4), a):
        for sg in permutations(range(a)):
            M = Ts[0][I[sg[0]]][I[0]]
            for k in range(1,a): M = M*Ts[k][I[sg[k]]][I[k]]
            acc += perm_sign(sg)*M
    return acc.applyfunc(sp.expand)
th = [d, 0]
t0=time.time()
vac = sp.zeros(1,16); vac[0, 0*4+3] = 1
A = [t(a, d-h, th) for a in range(4)]
Y = t(1, 0, th)
print('built', time.time()-t0, flush=True)
rows = []
for a in range(4):
    v = vac*A[a]
    for k in range(4):
        rows.append(v.applyfunc(sp.expand)); v = v*Y
M = sp.Matrix.vstack(*rows)
detM = sp.factor(M.det(method='berkowitz'))
print('det M(d) =', detM, time.time()-t0)
import pickle; pickle.dump({'M':M,'A':A,'Y':Y}, open('L2sym.pkl','wb'))
