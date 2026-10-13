"""Exact (Fraction) python cross-check of the Wa-sov-expansion-L4.wb results (written while the shared Wolfram kernel was busy).
Uses sov_l4.py machinery. Run from this folder. Caches pickles in ./cache."""
import sys, os, pickle, itertools, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'session-2026-10-12', 'sovq')))
from sov_l4 import *   # Fr, Chain, sov_words, ket_words, gram, measure, minv, Gc, vandermonde, WORDS, Z0, TH0, H0, product_vec, unit

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
os.makedirs(CACHE, exist_ok=True)


def cached(name, fn):
    f = os.path.join(CACHE, name + '.pkl')
    if os.path.exists(f):
        return pickle.load(open(f, 'rb'))
    r = fn()
    pickle.dump(r, open(f, 'wb'))
    return r


def charge_blocks():
    b = {}
    for i, x in enumerate(WORDS):
        b.setdefault(sum(x), []).append(i)
    return b


def comp_charge(c):   # charge of a product-basis component (0-based indices c = 4 digits)
    i1, i2, i3, i4 = c
    return i1 + (3 - i2) + i3 + (3 - i4)


COMPS = list(itertools.product(range(4), repeat=4))


def build(verbose=True):
    t0 = time.time()
    ch = Chain()
    W = cached('W', lambda: sov_words(ch))
    if verbose: print('W', time.time() - t0, flush=True)
    # inverse of W by blocks: block d: rows = words with |x| = d, cols = components with charge d
    def invW():
        Y = [[Fr(0)] * 256 for _ in range(256)]
        for d, xs in charge_blocks().items():
            cs = [j for j, c in enumerate(COMPS) if comp_charge(c) == d]
            B = [[W[i][j] for j in cs] for i in xs]
            Bi = minv(B)    # rows: comps, cols: words
            for a, j in enumerate(cs):
                for b, i in enumerate(xs):
                    Y[j][i] = Bi[a][b]
        return Y
    Y = cached('Y', invW)   # Y[comp][word]; W.Y = 1
    if verbose: print('Y', time.time() - t0, flush=True)
    Yk = cached('Yk', lambda: ket_words(ch))
    if verbose: print('Yk', time.time() - t0, flush=True)
    G = cached('Gram', lambda: gram(W, Yk))
    M = cached('M', lambda: measure(G))
    if verbose: print('M', time.time() - t0, flush=True)
    return ch, W, Y, Yk, G, M


def mat_rank(rows):
    """exact rank of a list of lists of Fractions"""
    A = [list(r) for r in rows]
    n = len(A); m = len(A[0]) if n else 0
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        for i in range(r + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / pv
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == n: break
    return r


def flat(c4, rows_idx, cols_idx):
    """flatten a 4-index tensor given as dict/array c4[(i,j,k,l)] into a matrix with given index groups"""
    import numpy as np
    out = []
    for ri in itertools.product(range(4), repeat=len(rows_idx)):
        row = []
        for ci in itertools.product(range(4), repeat=len(cols_idx)):
            idx = [0] * 4
            for p, v in zip(rows_idx, ri): idx[p] = v
            for p, v in zip(cols_idx, ci): idx[p] = v
            row.append(c4[tuple(idx)])
        out.append(row)
    return out


def tensor_ranks(vec):
    c4 = {w: vec[i] for i, w in enumerate(WORDS)}
    sites = [mat_rank(flat(c4, [s], [t for t in range(4) if t != s])) for s in range(4)]
    return {'nnz': sum(1 for x in vec if x != 0), 'site': sites,
            'pair12|34': mat_rank(flat(c4, [0, 1], [2, 3])), '13|24': mat_rank(flat(c4, [0, 2], [1, 3]))}
