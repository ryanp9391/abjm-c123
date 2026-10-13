import sys, time, pickle
from fractions import Fraction
import numpy as np
from modlib import *
from pipeline import build_positive_words, factor_list
from words import factor_list

def Wex_pair(a):
    return [Fraction(1), Fraction(0), Fraction(0), Fraction(a)], [Fraction(-a), Fraction(0), Fraction(0), Fraction(1)]

def kron_frac(vs):
    out = [Fraction(1)]
    for v in vs:
        out = [x * y for x in out for y in v]
    return out

def vec_mod(vs, p):
    x = kron_frac(vs)
    return np.array([fm(c, p) for c in x], dtype=np.float64)

def mod_targets(L, a, z, p):
    """exact targets as modular vectors (dict name -> vec)"""
    a = Fraction(a)
    v, u = Wex_pair(a)
    rest = []
    for _ in range(L // 2 - 1):
        vv, uu = Wex_pair(a); rest += [vv, uu]
    T = {}
    for a1 in (Fraction(-2, 3), Fraction(0), Fraction(2, 3), Fraction(5, 4)):
        vv, uu = Wex_pair(a1)
        T[f"T1 a'={a1}"] = vec_mod([vv, uu] + rest, p)
    for k in (Fraction(1, 3), Fraction(-2, 5), Fraction(5, 3)):
        vD = [v[i] * (1 - k * z[i]) for i in range(4)]
        uD = [u[i] / (1 - k * z[i]) for i in range(4)]
        T[f"T2 k={k}"] = vec_mod([vD, uD] + rest, p)
    for k in (Fraction(1, 3), Fraction(-2, 5), Fraction(5, 3)):
        vD = [v[i] * (1 - k / z[i]) for i in range(4)]
        T[f"T3a k={k}"] = vec_mod([vD, u] + rest, p)
    T["T3b vL^-1 x u(1+L)"] = vec_mod([[v[i] / z[i] for i in range(4)], [u[i] * (1 + z[i]) for i in range(4)]] + rest, p)
    T["T3b vL^-1 x u(2-L)"] = vec_mod([[v[i] / z[i] for i in range(4)], [u[i] * (2 - z[i]) for i in range(4)]] + rest, p)
    return T, vec_mod(rest, p) if rest else np.array([1.0])

def localised_dim_mod(basis, w0, p):
    dimv = basis.dim
    R = dimv // 16
    w0o = w0.astype(np.int64).astype(object)
    ww = int(np.dot(w0o, w0o) % p)
    winv = pow(ww, -1, p)
    # P = I - w0 w0^T / ww  (R x R) mod p
    P = (np.eye(R, dtype=object) - np.outer(w0o, w0o) * winv) % p
    P = P.astype(np.float64)
    rows = basis.rows.reshape(-1, 16, R)
    # rows x (16,R) @ P : products < 2^44, inner R up to 256 -> block
    out = np.zeros_like(rows)
    for b0 in range(0, R, 128):
        b1 = min(R, b0 + 128)
        out = np.mod(out + np.mod(rows[:, :, b0:b1] @ P[b0:b1, :], p), p)
    flat = out.reshape(basis.rank, -1)
    sub = ModBasis(flat.shape[1], p)
    sub.add(flat)
    return basis.rank - sub.rank

if __name__ == '__main__':
    L = int(sys.argv[1]); thp = int(sys.argv[2]); pi = int(sys.argv[3]); classes = tuple(sys.argv[4].split(','))
    p = PRIMES[pi]
    th = TH6 if thp == 1 else TH6b
    mc = ModChain(L, th, H, Z4, p)
    a0 = Fraction(1, 3)
    v, u = Wex_pair(a0)
    psi0 = vec_mod([v, u] * (L // 2), p)
    F = factor_list(mc)
    print('factors', len(F), flush=True)
    t = time.time()
    ws, B1, B2 = build_positive_words(mc, psi0, F, classes, verbose=False)
    W = ws.matrix()
    print('modular words built', W.shape, time.time() - t, flush=True)
    basis = ModBasis(mc.dim, p)
    t = time.time()
    # add in chunks, record for each class boundary the rank
    N = W.shape[1]
    lab = ws.labels
    nf = [len(l) for l in lab]
    cut = {}
    order = np.argsort(nf, kind='stable')
    Wn = W[:, order]
    nfs = np.array(nf)[order]
    res = {}
    for k in sorted(set(nf)):
        idx = np.nonzero(nfs == k)[0]
        for c0 in range(0, len(idx), 400):
            basis.add(Wn[:, idx[c0:c0 + 400]].T)
        res[k] = basis.rank
        print('words with <=', k, 'factors: rank', basis.rank, 'time', time.time() - t, flush=True)
        if k == 1 or k == 2 or True:
            T, w0 = mod_targets(L, a0, Z4, p)
            print('  loc dim', localised_dim_mod(basis, w0, p), ' membership:', {kk: basis.contains(vv) for kk, vv in T.items()}, flush=True)
    pickle.dump(dict(rank=res), open(f'modrun_L{L}_th{thp}_p{pi}.pkl', 'wb'))
