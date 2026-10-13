import numpy as np, itertools, time
from fractions import Fraction
from lib import *

MS = list(range(-4, 5))

def factor_list(ch, sites=(0, 1), ms=MS):
    F = []
    for s in sites:
        for a in (1, 2, 3):
            for m in ms:
                u = ch.th[s] + m * ch.h
                if ch.gprod(a, u) == 0:
                    continue
                F.append(dict(a=a, s=s, m=m, u=u, name=f"t{a}(th{s+1}{m:+d}h)"))
    return F

def apply_f(ch, f, Psi):
    return ch.apply(f['a'], f['u'], Psi)

def dense_t(ch, which, u, chunk=256):
    n = ch.dim
    M = np.zeros((n, n))
    for c0 in range(0, n, chunk):
        c1 = min(n, c0 + chunk)
        E = np.zeros((n, c1 - c0)); 
        for j in range(c0, c1): E[j, j - c0] = 1.0
        # rows of M: psi=e_j -> psi.t = row j
        M[c0:c1, :] = ch.apply(which, u, E).T
    return M

def colnorm(M):
    nr = np.linalg.norm(M, axis=0)
    nr[nr == 0] = 1.0
    return M / nr, nr
