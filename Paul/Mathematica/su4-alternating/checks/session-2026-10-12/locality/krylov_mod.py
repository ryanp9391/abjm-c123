import sys, time, pickle, random
from fractions import Fraction
import numpy as np
from modlib import *
from modrun import mod_targets, localised_dim_mod, Wex_pair, vec_mod

def run(L, thp, pi, seed=1, ustars=None, verbose=True):
    p = PRIMES[pi]
    th = TH6 if thp == 1 else TH6b
    mc = ModChain(L, th, H, Z4, p)
    a0 = Fraction(1, 3)
    v, u = Wex_pair(a0)
    psi0 = vec_mod([v, u] * (L // 2), p)
    rnd = random.Random(seed)
    # generic generator M = t1(u1) + c2 t2(u2) + c3 t3(u3)
    u1, u2, u3 = [Fraction(rnd.randint(-97, 97), rnd.randint(11, 53)) for _ in range(3)] if ustars is None else ustars
    c2, c3 = rnd.randint(2, p - 2), rnd.randint(2, p - 2)
    def M(x):
        y1 = mc.apply(1, u1, x)
        y2 = mc.apply(2, u2, x)
        y3 = mc.apply(3, u3, x)
        return np.mod(y1 + np.mod(c2 * y2, p) + np.mod(c3 * y3, p), p)
    basis = ModBasis(mc.dim, p)
    x = psi0.copy()
    t0 = time.time()
    K = 0
    while True:
        fl = basis.add(x[None, :])
        if not fl[0]:
            break
        K += 1
        if verbose and K % 50 == 0:
            print('krylov step', K, 'time', time.time() - t0, flush=True)
        x = M(x)
    T, w0 = mod_targets(L, a0, Z4, p)
    mem = {k: basis.contains(vv) for k, vv in T.items()}
    ld = localised_dim_mod(basis, w0, p)
    return K, ld, mem

if __name__ == '__main__':
    L = int(sys.argv[1]); thp = int(sys.argv[2]); pi = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    K, ld, mem = run(L, thp, pi, seed)
    print('L', L, 'theta pt', thp, 'prime', PRIMES[pi], 'seed', seed)
    print('cyclic (Krylov) dimension', K, ' localised dim', ld)
    for k, v in mem.items(): print('  ', k, v)
