# ABJM point, binary-exact floats: span of ALL individually localising words (products up to length 4 of
# t_k(theta_F + m h), k=1..3, m=-2..5, regularised) on generic BPS references, random twist
import numpy as np, itertools, sys, time
from pgbase import *
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 3
rng = np.random.default_rng(seed)
h = 0.5; thF = -0.625; thA = thF + 2*h; th3, th4 = 0.8125+1.0, 0.8125
mode = sys.argv[2] if len(sys.argv) > 2 else 'rand'
G = rng.normal(size=(4,4)) if mode == 'rand' else np.diag([1.5,1.25,0.625,0.5625])
ch = Chain(['A','F','A','F'], [thA, thF, th3, th4], h, G.astype(object))
if mode == 'Wa':
    a = 0.375; v = np.array([1,0,0,a]); u = np.array([-a,0,0,1.0])
else:
    v = rng.normal(size=4); u = rng.normal(size=4); u = u - (v@u)/(v@v)*v
ref = prod_state([v,u,v,u], dtype=float); rest = np.kron(v,u)
def nonloc(out):
    M = out.reshape(16,16); n = np.linalg.norm(M)
    if n < 1e-10: return -1, None
    ps = M@rest/(rest@rest)
    return np.linalg.norm(M - np.outer(ps, rest))/n, ps
ops = [(k, m) for k in (1,2,3) for m in range(-2, 6)]
cache = {(): ref}
loc = {(): rest*0 + 0}
loc[()] = np.kron(v,u)
t0 = time.time()
for L in (1,2,3,4):
    for w in itertools.combinations_with_replacement(range(len(ops)), L):
        prev = cache.get(w[:-1])
        if prev is None: continue
        k, m = ops[w[-1]]
        out = ch.apply_t(k, thF + m*h, prev, reg=True)
        out = out/ max(np.abs(out).max(), 1e-300)
        if L < 4: cache[w] = out
        r, ps = nonloc(out)
        if 0 <= r < 1e-9: loc[w] = ps/np.linalg.norm(ps)
    print('length', L, 'done, #localising so far', len(loc), 'time %.0f' % (time.time()-t0)); sys.stdout.flush()
B = np.array(list(loc.values())).T
sv = np.linalg.svd(B, compute_uv=False)
print('rank of localised span:', np.sum(sv > 1e-8*sv[0]), np.round(sv[:8], 10))
for w in list(loc)[:60]: print([ops[i] for i in w])
