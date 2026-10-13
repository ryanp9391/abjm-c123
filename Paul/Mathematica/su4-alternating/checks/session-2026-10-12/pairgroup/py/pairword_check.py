# check: localised pair-1 state of each localising word (4 and 6 sites) == lambda * (same word on the 1-pair chain)
import numpy as np, itertools
from pgbase import *
h = Fr(1,2); thF = Fr(-3,5); thA = thF + 2*h
th4 = [thA, thF, Fr(-3,8), Fr(-11,8)]
th6 = th4 + [Fr(5,13)+2*h, Fr(5,13)]
rng = np.random.default_rng(2)
G = np.array([[Fr(int(rng.integers(-5,6)), int(rng.integers(1,4))) for j in range(4)] for i in range(4)], dtype=object)
v = [Fr(1), Fr(2), Fr(-1), Fr(3)]; u = [Fr(2), Fr(1,2), Fr(1), Fr(0)]
uv = sum(x*y for x,y in zip(v,u)); u = [u[i] - uv/sum(x*x for x in v)*v[i] for i in range(4)]
assert sum(x*y for x,y in zip(v,u)) == 0
words = [[(1,0)], [(2,1)], [(3,1)], [(1,0),(1,1)], [(1,0),(3,1)], [(1,0),(3,3)], [(1,2),(3,1)], [(3,1),(3,2)],
         [(1,0),(1,0),(3,3)], [(1,0),(1,2),(3,1)], [(1,0),(2,1),(3,3)], [(1,0),(2,2),(3,1)], [(1,0),(3,1),(3,3)], [(1,2),(2,1),(3,1)], [(1,2),(3,1),(3,1)]]
def apply_word(ch, w, psi):
    for k, m in w: psi = ch.apply_t(k, thF + m*h, psi, reg=True)
    return psi
chP = Chain(['A','F'], [thA, thF], h, G)
ref2 = prod_state([v,u]); ref2 = np.vectorize(Fr, otypes=[object])(ref2)
for N, th in ((4, th4),):
    ch = Chain(['A','F']*(N//2), th, h, G)
    ref = prod_state([v,u]*(N//2)); ref = np.vectorize(Fr, otypes=[object])(ref)
    rest = prod_state([v,u]*(N//2-1)).reshape(-1)
    rest = np.vectorize(Fr, otypes=[object])(rest)
    rr = sum(x*x for x in rest)
    for w in words:
        out = apply_word(ch, w, ref).reshape(16, -1)
        pairst = out.dot(rest) / rr
        nonloc = out - np.multiply.outer(pairst, rest)
        nl = max(abs(float(x)) for x in nonloc.reshape(-1))
        pw = apply_word(chP, w, ref2).reshape(-1)
        # ratio
        idx = [i for i in range(16) if pw[i] != 0]
        if not idx:
            print(N, w, 'pair word zero; full nonloc', nl, 'pair state zero?', all(x == 0 for x in pairst)); continue
        lam = pairst[idx[0]] / pw[idx[0]]
        diff = max(abs(float(pairst[i] - lam*pw[i])) for i in range(16))
        print(N, w, 'nonloc', nl, ' pair-word match diff', diff, ' lambda', lam)
