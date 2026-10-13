# Exact full-chain verification (Fractions) of the pair word K_kappa at L = 2, 3 pairs:
#   <W_a^L| K_kappa = (1 - kappa z1)^2 <W_{a'} (x) W_a^{L-1}|,  ABJM point, z2 = z3 = y
import numpy as np, sys, time
from pgbase import *
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
h = Fr(1,2)
thF = [Fr(-3,5), Fr(-11,8), Fr(-7,9)][:L]
th = []
for f in thF: th += [f + 2*h, f]
y = Fr(5,4); z = [Fr(3,2), y, y, Fr(4,7)]
ch = Chain(['A','F']*L, th, h, diagm(z))
def lam(k, u0):
    r = Fr(1)
    for s in range(2, 2*L):
        t = th[s]
        if s % 2 == 0:  # A
            for j in range(k): r *= (u0 - t - j*h)
        else:
            r *= (u0 - t - k*h)
            for j in range(k-1): r *= (u0 - t - j*h)
    return r
u1 = th[1]; u2 = th[1] + h
l1, l2, l3 = lam(1, u1), lam(2, u2), lam(3, u2)
print('L =', L, ' lambda1, lambda2, lambda3 =', l1, l2, l3)
def Wp(a): return [[Fr(1),Fr(0),Fr(0),a],[-a,Fr(0),Fr(0),Fr(1)]]
def run(a, kaps, ref=None, target_pair=None):
    if ref is None: ref = prod_state(Wp(a)*L)
    t0 = time.time()
    x1 = ch.apply_t(1, u1, ref)
    x2 = ch.apply_t(2, u2, ref)
    x3 = ch.apply_t(3, u2, ref, reg=True)
    x31 = ch.apply_t(1, u1, x3)
    print('   (4 transfer-matrix applications: %.1f s)' % (time.time()-t0))
    for kap in kaps:
        s = 1 - kap*(z[0] + z[3])
        K = s*s*ref + kap*s*(x1/(h**2*l1) + x2/(2*h**4*y*l2) + x3/(6*h**5*y**2*l3)) + kap**2*x31/(24*h**7*y**2*l1*l3)
        ap = a*(1-kap*z[3])/(1-kap*z[0])
        tgt = (1-kap*z[0])**2 * prod_state(Wp(ap) + Wp(a)*(L-1))
        diff = K - tgt
        nz = sum(1 for x in diff.reshape(-1) if x != 0)
        print('   a = %s, kappa = %s: a\' = %s, nonzero entries of difference: %d' % (a, kap, ap, nz))
run(Fr(1,3), [Fr(3,10), Fr(-7,5), 1/z[3]])
run(Fr(-1), [Fr(2,9), 1/z[3]])
