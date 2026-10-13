# Exact negative controls at 4 sites (Fractions, sympy rank):
#  (1) z2 != z3, ABJM point, W_a^2: rank(localised span) vs rank(+target)
#  (2) z2 == z3, generic theta (non-ABJM): same
#  (3) z2 == z3, ABJM point, generic BPS reference (v (x) u)^2: span vs group-rotated pair, and the K word itself
import numpy as np, sympy as sp
from pgbase import *
h = Fr(1,2)
def setup(thF1, d, z, thF2=Fr(-11,8)):
    th = [thF1 + d, thF1, thF2 + 2*h, thF2]
    return Chain(['A','F','A','F'], th, h, diagm(z)), th
def pairstates(ch, th, ref, rest, words):
    rr = rest.dot(rest); cols = []
    for w in words:
        psi = ref
        for (k, m) in w: psi = ch.apply_t(k, th[1] + m*h, psi, reg=True)
        M = psi.reshape(16, 16); ps = M.dot(rest)/rr
        nl = M - np.multiply.outer(ps, rest)
        cols.append((ps, all(x == 0 for x in nl.reshape(-1))))
    return cols
W5 = [[], [(1,0)], [(2,1)], [(3,1)], [(1,0),(3,1)]]
Wg = [[], [(1,0)], [(3,-1)], [(1,0),(3,-1)]]       # generic theta: t1(th2), t3(th1-h) = t3(th2 + d - h)
def Wp(a): return [[Fr(1),Fr(0),Fr(0),a],[-a,Fr(0),Fr(0),Fr(1)]]
def rk(cols, tgt=None):
    B = sp.Matrix(16, len(cols), lambda i, j: cols[j][i])
    if tgt is None: return B.rank()
    return B.rank(), B.row_join(sp.Matrix(list(tgt))).rank()
a = Fr(1,3); kap = Fr(3,10)
def target(v, u, z, kap):
    Dm = [1 - kap*zi for zi in z]
    vD = [v[i]*Dm[i] for i in range(4)]; uD = [u[i]/Dm[i] for i in range(4)]
    return prod_state([vD, uD]).reshape(-1)
# (1)
z = [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
ch, th = setup(Fr(-3,5), 2*h, z)
ref = prod_state(Wp(a)*2); rest = prod_state(Wp(a)).reshape(-1)
cols = pairstates(ch, th, ref, rest, W5)
print('(1) ABJM, z2!=z3, W_a^2: all localise:', all(c[1] for c in cols), ' rank span, rank +target:', rk([c[0] for c in cols], target(Wp(a)[0], Wp(a)[1], z, kap)))
# (2)
z = [Fr(3,2), Fr(5,4), Fr(5,4), Fr(4,7)]
ch, th = setup(Fr(-3,5), Fr(2,7) + Fr(3,5), z)
words_g = [[], [(1,0)], [(3, None)], [(1,0),(3,None)]]
cols = []
rr = rest.dot(rest)
for w in words_g:
    psi = ref
    for (k, m) in w:
        u0 = th[1] if k == 1 else th[0] - h
        psi = ch.apply_t(k, u0, psi, reg=True)
    M = psi.reshape(16,16); ps = M.dot(rest)/rr
    cols.append((ps, all(x == 0 for x in (M - np.multiply.outer(ps, rest)).reshape(-1))))
print('(2) generic theta, z2=z3, W_a^2 words {1,t1(th2),t3(th1-h),prod}: all localise:', all(c[1] for c in cols), ' rank span, rank +target:', rk([c[0] for c in cols], target(Wp(a)[0], Wp(a)[1], z, kap)))
# also try adding t2(th2+h) at generic theta (does not localise) -> report nonlocality
psi = ch.apply_t(2, th[1] + h, ref); M = psi.reshape(16,16); ps = M.dot(rest)/rr
print('    t2(th2+h) at generic theta localises?', all(x == 0 for x in (M - np.multiply.outer(ps, rest)).reshape(-1)))
# (3)
ch, th = setup(Fr(-3,5), 2*h, z)
v = [Fr(1), Fr(2), Fr(-1), Fr(3)]; u0 = [Fr(2), Fr(1,2), Fr(1), Fr(1)]
vu = sum(x*y for x, y in zip(v, u0)); vv = sum(x*x for x in v)
u = [u0[i] - vu/vv*v[i] for i in range(4)]
ref3 = prod_state([v, u, v, u]); rest3 = prod_state([v, u]).reshape(-1)
cols = pairstates(ch, th, ref3, rest3, W5)
print('(3) ABJM, z2=z3, generic BPS: all localise:', all(c[1] for c in cols), ' rank span, rank +target:', rk([c[0] for c in cols], target(v, u, z, kap)))
# the K word on the generic BPS reference: product state? compare with target
y = z[1]
def lam(k, u0):
    r = Fr(1)
    for s in (2, 3):
        t = th[s]
        fac = [u0 - t - j*h for j in range(k)] if s == 2 else [u0 - t - k*h] + [u0 - t - j*h for j in range(k-1)]
        for f in fac:
            if f != 0: r *= f
    return r
l1, l2, l3 = lam(1, th[1]), lam(2, th[1]+h), lam(3, th[1]+h)
s = 1 - kap*(z[0] + z[3])
Kst = s*s*cols[0][0] + kap*s*(cols[1][0]/(h**2*l1) + cols[2][0]/(2*h**4*y*l2) + cols[3][0]/(6*h**5*y**2*l3)) + kap**2*cols[4][0]/(24*h**7*y**2*l1*l3)
Km = sp.Matrix(4, 4, lambda i, j: Kst[4*i+j])
tg = target(v, u, z, kap); Tm = sp.Matrix(4, 4, lambda i, j: tg[4*i+j])*(1-kap*z[0])*(1-kap*z[3])
print('    K word on generic BPS pair: rank of 4x4 pair matrix =', Km.rank(), ' (product state iff 1);  ||K - d1 d4 rho(D)||/||.|| =',
      float(sp.sqrt(sum((Km - Tm).applyfunc(lambda x: x**2))) / sp.sqrt(sum(Tm.applyfunc(lambda x: x**2)))))
sing = sum(Km[i,i] for i in range(4))
print('    trace of pair matrix (singlet admixture, BPS iff 0):', sing, '=', float(sing))
