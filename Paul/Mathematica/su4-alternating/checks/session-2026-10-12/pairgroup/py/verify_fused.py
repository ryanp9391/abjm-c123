# verify fused single-pass t_k against the quantum-minor definition (numeric float, N=2 and N=4)
import numpy as np, itertools
from pgbase import *
rng = np.random.default_rng(1)
def full_T(sites, th, h, G, u):
    """monodromy T^G(u) as array [i,j, M, M] (operator matrices acting on bras from the right)"""
    N = len(sites); dim = 4**N
    T = None
    for s in range(N):
        L = lax_k(sites[s], 1, u, th[s], h).astype(float)  # [i,j,p,q]
        if T is None: T = L
        else:
            # T[i,k] (x) L[k,j]
            T = np.einsum('ikab,kjpq->ijapbq', T, L).reshape(4, 4, T.shape[2]*4, T.shape[3]*4)
    TG = np.einsum('ikab,kj->ijab', T, G)
    return TG
def tk_def(sites, th, h, G, u, k):
    N = len(sites); dim = 4**N
    Ts = [full_T(sites, th, h, G, u - m*h) for m in range(k)]
    tot = np.zeros((dim, dim))
    for I in itertools.combinations(range(4), k):
        for p in itertools.permutations(range(k)):
            s = perm_sign(p)
            M = np.eye(dim)
            for m in range(k):
                M = M @ Ts[m][I[p[m]], I[m]]
            tot += s * M
    return tot
for sites, th in ((['A','F'], [0.2857, -0.6]), (['A','F','A','F'], [2/7, -3/5, 8/11, -11/8]), (['F','A'], [0.3,-0.45])):
    h = 0.5
    G = rng.normal(size=(4,4))
    ch = Chain(sites, th, h, G.astype(object))
    N = len(sites)
    for k in (1,2,3,4):
        u = 0.37
        Mdef = tk_def(sites, th, h, G, u, k)
        # single pass on random bras
        err = 0
        for trial in range(3):
            psi = rng.normal(size=(4,)*N)
            a = ch.apply_t(k, u, psi.astype(float)).reshape(-1)
            b = psi.reshape(-1) @ Mdef
            err = max(err, np.abs(a-b).max()/np.abs(b).max())
        print(sites, 'k=',k, 'rel err', err)
