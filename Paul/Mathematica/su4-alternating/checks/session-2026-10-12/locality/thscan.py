import random, time
from fractions import Fraction
import numpy as np
from modlib import *
from modrun import Wex_pair, vec_mod, mod_targets
from cyc import ModBasis
p = PRIMES[0]
def kry(L, th, h, z, a0, seed=1):
    mc = ModChain(L, th, h, z, p)
    v, u = Wex_pair(a0)
    psi0 = vec_mod([v, u] * (L // 2), p)
    rnd = random.Random(seed)
    u1, u2, u3 = [Fraction(rnd.randint(-97, 97), rnd.randint(11, 53)) for _ in range(3)]
    c2, c3 = rnd.randint(2, p - 2), rnd.randint(2, p - 2)
    M = lambda x: np.mod(mc.apply(1, u1, x) + np.mod(c2 * mc.apply(2, u2, x), p) + np.mod(c3 * mc.apply(3, u3, x), p), p)
    B = ModBasis(mc.dim, p); x = psi0.copy()
    while B.add(x[None, :])[0]: x = M(x)
    return B.rank
rnd = random.Random(5)
for trial in range(6):
    th = [Fraction(rnd.randint(-60, 60), rnd.randint(7, 40)) for _ in range(4)]
    print('random theta', [str(x) for x in th], 'L=2 cyclic dim', kry(4, th, H, Z4, Fraction(1, 3)))
print('std', kry(4, TH6, H, Z4, Fraction(1, 3)), 'theta2', kry(4, TH6b, H, Z4, Fraction(1, 3)))
# dependence on a and z
for a0 in (Fraction(1, 3), Fraction(2, 5), Fraction(-3, 7)):
    print('a=', a0, kry(4, TH6, H, Z4, a0))
print('z2', kry(4, TH6, H, [Fraction(2), Fraction(3, 5), Fraction(7, 3), Fraction(5, 4)], Fraction(1, 3)))
