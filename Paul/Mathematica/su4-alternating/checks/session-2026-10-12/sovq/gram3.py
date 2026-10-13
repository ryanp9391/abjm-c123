import pickle, itertools
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from sov_l4 import WORDS, minv
d = pickle.load(open('gram_exact.pkl', 'rb')); G = d['G']
idx = {w: i for i, w in enumerate(WORDS)}
bar = lambda y: tuple(3 - v for v in y)
# Gbar[x][x'] = <x| ybar=x'>  where ket label y = 3 - x'
Gb = {(x, xp): G[idx[x], idx[bar(xp)]] for x in WORDS for xp in WORDS if G[idx[x], idx[bar(xp)]] != 0}
print('diag nonzero:', all((x, x) in Gb for x in WORDS))
# triangularity: directed graph x -> x' for off-diagonal nonzeros; acyclic?
adj = defaultdict(list)
for (x, xp) in Gb:
    if x != xp: adj[x].append(xp)
state = {}
def dfs(u):
    state[u] = 1
    for v in adj[u]:
        if state.get(v) == 1: return True
        if state.get(v) is None and dfs(v): return True
    state[u] = 2
    return False
import sys; sys.setrecursionlimit(10000)
cyc = any(dfs(u) for u in WORDS if state.get(u) is None)
print('off-diagonal graph has a cycle (not triangularisable):', cyc)
# symmetric pattern?
sym = sum(1 for (x, xp) in Gb if (xp, x) in Gb)
print('# nonzeros', len(Gb), '# with transposed partner nonzero', sym)
# per-pair (A1F2),(A3F4) factorisation test of the pattern: Gb[(x,xp)] != 0 iff P12(x12,xp12) and P34(...) ?
p12 = set((x[:2], xp[:2]) for (x, xp) in Gb); p34 = set((x[2:], xp[2:]) for (x, xp) in Gb)
fact = sum(1 for x in WORDS for xp in WORDS if (x[:2], xp[:2]) in p12 and (x[2:], xp[2:]) in p34 and sum(x) == sum(xp))
print('pattern from pair projections with charge conservation:', fact, 'vs actual', len(Gb))
# blocks by charge; exact inverse = measure
blocks = defaultdict(list)
for x in WORDS: blocks[sum(x)].append(x)
M = {}
for S, xs in sorted(blocks.items()):
    B = [[Gb.get((x, xp), Fr(0)) for xp in xs] for x in xs]
    Bi = minv(B)
    for i, x in enumerate(xs):
        for j, xp in enumerate(xs):
            if Bi[i][j] != 0: M[(x, xp)] = Bi[i][j]   # key (xp_label, x_label): row i of B^-1 is xp, column j is x
    nnzB = sum(1 for r in B for v in r if v != 0); nnzI = sum(1 for r in Bi for v in r if v != 0)
    print('charge', S, 'block', len(xs), 'nnz G', nnzB, 'nnz G^-1', nnzI)
pickle.dump({'Gb': Gb, 'M': M}, open('measure_exact.pkl', 'wb'))
