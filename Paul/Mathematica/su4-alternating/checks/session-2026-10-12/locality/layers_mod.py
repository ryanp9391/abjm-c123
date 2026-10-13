import sys, time, pickle
from fractions import Fraction
import numpy as np
from modlib import *
from modrun import mod_targets, localised_dim_mod, Wex_pair, vec_mod
from words import factor_list

def layered(L, thp, pi, base_vecs_fn=None, maxlayers=8, factors=None, tag='', verbose=True):
    p = PRIMES[pi]
    th = TH6 if thp == 1 else TH6b
    mc = ModChain(L, th, H, Z4, p)
    a0 = Fraction(1, 3)
    v, u = Wex_pair(a0)
    psi0 = vec_mod([v, u] * (L // 2), p)
    F = factor_list(mc) if factors is None else factors(mc)
    T, w0 = mod_targets(L, a0, Z4, p)
    basis = ModBasis(mc.dim, p)
    basis.add(psi0[None, :])
    new = psi0[None, :]
    out = []
    t0 = time.time()
    def rec(l):
        ld = localised_dim_mod(basis, w0, p)
        mem = {k: basis.contains(vv) for k, vv in T.items()}
        out.append((l, basis.rank, ld, mem))
        if verbose:
            print(f'layer {l}: rank {basis.rank}, loc {ld}, T1 {[mem[k] for k in mem if k.startswith("T1")]}, T2 {[mem[k] for k in mem if k.startswith("T2")]}, T3 {[mem[k] for k in mem if k.startswith("T3")]}  t={time.time()-t0:.0f}s', flush=True)
    rec(0)
    for l in range(1, maxlayers + 1):
        cand = []
        for f in F:
            cand.append(mc.apply(f['a'], f['u'], new.T).T)
        cand = np.concatenate(cand, axis=0)
        # add in chunks, collect independent rows
        added = []
        for c0 in range(0, cand.shape[0], 400):
            ch = np.mod(cand[c0:c0 + 400], p)
            fl = basis.add(ch)
            added.append(ch[np.array(fl)])
        new = np.concatenate(added, axis=0) if added else np.zeros((0, mc.dim))
        rec(l)
        if new.shape[0] == 0:
            break
    return out

if __name__ == '__main__':
    L = int(sys.argv[1]); thp = int(sys.argv[2]); pi = int(sys.argv[3])
    out = layered(L, thp, pi)
    pickle.dump(out, open(f'layers_L{L}_th{thp}_p{pi}.pkl', 'wb'))
