# single 6-site: fused wedge^2-aux Lax from quantum minors of L(u)L(u-h); eigen-structure vs u
import numpy as np, itertools, sympy as sp
from pgbase import subsets, DK, perm_sign
u = sp.Symbol('u'); h = sp.Rational(1,2)
S2 = subsets(2)
D2 = {(a,b): sp.Matrix(DK[2][(a,b)].tolist()) for a in range(4) for b in range(4)}
def L(uu):
    return [[(uu if i == j else 0)*sp.eye(6) + h*D2[(i,j)] for j in range(4)] for i in range(4)]
L0 = L(u); L1 = L(u - h)
F = sp.zeros(36, 36)
for x, I in enumerate(S2):
    for y, J in enumerate(S2):
        # qmin_{I,J} = sum_sigma sgn T_{I_s1, J1}(u) T_{I_s2, J2}(u-h)
        M = sp.zeros(6,6)
        for p in itertools.permutations(range(2)):
            M += perm_sign(p) * L0[I[p[0]]][J[0]] * L1[I[p[1]]][J[1]]
        F[x*6:(x+1)*6, y*6:(y+1)*6] = M
F = F.applyfunc(sp.expand)
# eigenvalues as functions of u (block structure: aux (x) site, both 6-dim)
ev = F.eigenvals()
for k, m in ev.items(): print(sp.factor(k), 'mult', m)
# permutation operator on 6(x)6
P = sp.zeros(36,36)
for x in range(6):
    for y in range(6): P[x*6+y, y*6+x] = 1
for val in [-2*h, -3*h/2, -h, -h/2, 0, h/2, h]:
    Fv = F.subs(u, val)
    # is Fv proportional to P ?
    c = None
    ok = True
    for i in range(36):
        for j in range(36):
            if P[i,j] != 0:
                if c is None: c = Fv[i,j]
    print('u=', val, ' F - c P zero?', (Fv - c*P).is_zero_matrix, ' c=', c, ' rank F', Fv.rank())
