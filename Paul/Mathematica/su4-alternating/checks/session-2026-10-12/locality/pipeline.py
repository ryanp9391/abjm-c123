import numpy as np, itertools, time, pickle, os
import scipy.linalg as sla
from words import *

def get_inverse_lu(ch, f):
    M = dense_t(ch, f['a'], f['u'])
    lu = sla.lu_factor(M)
    d = np.abs(np.diag(lu[0]))
    return lu, d.min() / d.max()

def apply_inv(lu, Psi):
    # Psi: (dim,) or (dim,B): returns psi . M^{-1} columnwise
    return sla.lu_solve(lu, Psi, trans=1)

class WordSet:
    def __init__(self):
        self.vecs = []   # list of arrays (dim, k)
        self.labels = []  # list of tuples of ((fi,e),...) per column
    def add(self, V, labs):
        self.vecs.append(V)
        self.labels.extend(labs)
    def matrix(self):
        return np.concatenate(self.vecs, axis=1)

def build_positive_words(ch, psi0, F, classes=('W1', 'W2', 'W3AFF', 'W3AAF', 'W4'), verbose=True):
    """positive-power words; returns WordSet (identity word included)."""
    ws = WordSet()
    n = len(F)
    ws.add(psi0[:, None], [()])
    B1 = np.stack([apply_f(ch, f, psi0) for f in F], axis=1)
    ws.add(B1, [((i, 1),) for i in range(n)])
    Aidx = [i for i in range(n) if F[i]['s'] == 0]
    Fidx = [i for i in range(n) if F[i]['s'] == 1]
    # level 2 : all f<=f'
    B2 = {}
    if any(c in classes for c in ('W2', 'W3AFF', 'W3AAF')):
        t = time.time()
        for i in range(n):
            R = apply_f(ch, F[i], B1[:, i:])   # f_i applied to B1[j], j>=i  -> word (j, i)
            for jj, j in enumerate(range(i, n)):
                B2[(i, j)] = R[:, jj]
        if verbose: print('level2 done', time.time() - t, flush=True)
    if 'W2' in classes:
        keys = sorted(B2)
        ws.add(np.stack([B2[k] for k in keys], axis=1), [((k[0], 1), (k[1], 1)) for k in keys])
    # level 3: AFF
    if 'W3AFF' in classes:
        t = time.time()
        keys = [(i, j) for (i, j) in B2 if (i in Fidx and j in Fidx)]  # F<=F pairs
        keys.sort()
        P = np.stack([B2[k] for k in keys], axis=1)
        labs = []; blocks = []
        for ia in Aidx:
            R = apply_f(ch, F[ia], P)
            blocks.append(R); labs += [((ia, 1), (k[0], 1), (k[1], 1)) for k in keys]
        ws.add(np.concatenate(blocks, axis=1), labs)
        if verbose: print('W3AFF done', time.time() - t, flush=True)
    if 'W3AAF' in classes:
        t = time.time()
        keys = [(i, j) for (i, j) in B2 if (i in Aidx and j in Aidx)]
        keys.sort()
        P = np.stack([B2[k] for k in keys], axis=1)
        labs = []; blocks = []
        for iF in Fidx:
            R = apply_f(ch, F[iF], P)
            blocks.append(R); labs += [((k[0], 1), (k[1], 1), (iF, 1)) for k in keys]
        ws.add(np.concatenate(blocks, axis=1), labs)
        if verbose: print('W3AAF done', time.time() - t, flush=True)
    return ws, B1, B2

def factor_index(F, a, s, m):
    for i, f in enumerate(F):
        if f['a'] == a and f['s'] == s and f['m'] == m:
            return i
    return None


def build_inverse_words(ch, F, inv_idx, base_vecs, base_labels, verbose=True):
    """words (f,-1) x base words.  inv_idx: indices into F.  Returns (vecs, labels, info)."""
    vecs = []; labs = []; info = []
    for i in inv_idx:
        t = time.time()
        lu, rc = get_inverse_lu(ch, F[i])
        info.append((F[i]['name'], rc))
        if rc < 1e-13:
            if verbose: print('skip (near singular)', F[i]['name'], rc)
            continue
        V = apply_inv(lu, base_vecs)
        vecs.append(V); labs += [((i, -1),) + l for l in base_labels]
        if verbose: print('inv', F[i]['name'], 'pivot ratio %.2e' % rc, 'time %.1f' % (time.time() - t), flush=True)
    return (np.concatenate(vecs, axis=1) if vecs else np.zeros((ch.dim, 0))), labs, info


def label_str(F, lab):
    if not lab: return '1'
    return ' '.join(F[i]['name'] + ('' if e == 1 else '^-1') for i, e in lab)
