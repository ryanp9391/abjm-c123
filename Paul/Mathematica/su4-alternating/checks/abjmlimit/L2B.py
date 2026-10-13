import sympy as sp, pickle
from itertools import permutations, combinations
from chain import perm_sign
exec(open('L2num2.py').read().split('results = []')[0])
from L2sym import E, rho, types
d0 = 2*h
def Tun(uv):
    L = [[[ (uv-[d0,0][s])*(1 if i==j else 0)*sp.eye(4) + h*rho(types[s],i,j) for j in range(4)] for i in range(4)] for s in range(2)]
    return [[sum((sp.kronecker_product(L[0][i][k], L[1][k][j]) for k in range(4)), sp.zeros(16,16)) for j in range(4)] for i in range(4)]
def qmin(rows, cols, uv):
    a = len(cols); Ts = [Tun(uv - k*h) for k in range(a)]; acc = sp.zeros(16,16)
    for sg in permutations(range(a)):
        M = Ts[0][rows[sg[0]]][cols[0]]
        for k in range(1,a): M = M*Ts[k][rows[sg[k]]][cols[k]]
        acc += perm_sign(sg)*M
    return acc
def Bop(uv):
    acc = sp.zeros(16,16)
    for k1 in (1,2,3):
        for K2 in combinations((1,2,3),2):
            acc += Tun(uv)[k1-1][0]*qmin(tuple(k-1 for k in K2), (0,k1), uv+h)*qmin((0,1,2),(0,K2[0],K2[1]), uv+2*h)
    return acc
B = Bop(sp.Rational(37,100))
R = pickle.load(open('L2reg0.pkl','rb'))['R']
# Jordan form: eigenvalues and block sizes
ev = B.eigenvals()
from collections import Counter
blocks = []
for lam, mult in ev.items():
    r1 = (B - lam*sp.eye(16)).rank(); r2 = ((B - lam*sp.eye(16))**2).rank()
    blocks.append((mult, 16-r1, 16-r2))
print('# distinct eigenvalues of B (d=2h, L=2):', len(ev), ' (alg mult, geom mult, dim ker^2):', Counter(blocks))
# R rows: eigen or Jordan partner?
keys = [(a,k) for a in range(4) for k in range(4)]
for i,(a,k) in enumerate(keys):
    v = R.row(i); vB = v*B
    if sp.Matrix.vstack(v, vB).rank() == 1: continue
    # find lam with vB - lam v in span of eigen-row of word (1,k+1)
    j = keys.index((1,k+1))
    lam = None
    for l in ev:
        rr = vB - l*v
        if sp.Matrix.vstack(R.row(j), rr).rank() == 1 and not rr.is_zero_matrix:
            c = [rr[t]/R[j,t] for t in range(16) if R[j,t] != 0][0]
            print('  reg row', (a,k), ': v B = beta v + (%s) <1,%d|' % (sp.nsimplify(c), k+1))
