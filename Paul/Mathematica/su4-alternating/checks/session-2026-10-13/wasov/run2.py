"""W_a components, twisted vacuum <0_w|, global Borel rotations: exact python cross-check of the notebook."""
from wa_sov_core import *
import numpy as np, random
ch, W, Y, Yk, G, M = build(verbose=False)
Yn = np.array(Y, dtype=object); Ykn = np.array(Yk, dtype=object)
V = vandermonde(Z0); Vi = minv(V)
vz = [[V[i][j] for i in range(4)] for j in range(4)]
ut = [Vi[j] for j in range(4)]

def prod_bra(vs): return product_vec(vs)
def components(A0, A1, F0, F1):
    comp = {n: np.array([Fr(0)] * 256, dtype=object) for n in range(5)}
    for flips in itertools.product([0, 1], repeat=4):
        vs = [(A1 if fl else A0) if s % 2 == 0 else (F1 if fl else F0) for s, fl in enumerate(flips)]
        comp[sum(flips)] = comp[sum(flips)] + prod_bra(vs)
    return comp

# ---- W_a components: A0 = v(z1), A1 = v(z4), F0 = u4, F1 = -u1
comp = components(vz[0], vz[3], ut[3], [-x for x in ut[0]])
cn = {n: comp[n].dot(Yn) for n in comp}
print('== <W_a| components (n = power of a):')
for n in range(5):
    print(n, tensor_ranks(list(cn[n])))
# m(r2,r4) matrix rank at random points, for each component and the generic bra a = 3/7
def mrank(c, r2, r4):
    c4 = {w: c[i] for i, w in enumerate(WORDS)}
    return mat_rank([[sum(c4[(a1, k2, a3, k4)] * r2 ** k2 * r4 ** k4 for k2 in range(4) for k4 in range(4)) for a3 in range(4)] for a1 in range(4)])
random.seed(1)
pts = [(Fr(random.randint(2, 40), random.randint(3, 17)), Fr(random.randint(-30, -2), random.randint(3, 19))) for _ in range(3)]
for n in range(5):
    print('m rank comp', n, [mrank(list(cn[n]), *p) for p in pts])
a0 = Fr(3, 7)
cgen = sum((a0 ** n * cn[n] for n in range(5)), np.array([Fr(0)] * 256, dtype=object))
print('generic a=3/7', tensor_ranks(list(cgen)), [mrank(list(cgen), *p) for p in pts])
pickle.dump({n: list(cn[n]) for n in cn}, open(os.path.join(CACHE, 'cn.pkl'), 'wb'))

# ---- dual side: <W_a^{(n)}|y> tensors
print('== <W^(n)|y> tensor ranks (site ranks all 1 <=> product over sites):')
for n in range(5):
    print(n, tensor_ranks(list(comp[n].dot(Ykn))))
bgen = sum((a0 ** n * comp[n] for n in range(5)), np.array([Fr(0)] * 256, dtype=object)).dot(Ykn)
print('generic a=3/7', tensor_ranks(list(bgen)))

# ---- twisted vacuum <0_w|
Wt = (Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9))
Vw = vandermonde(Wt); Viw = minv(Vw)
vw = [Vw[i][0] for i in range(4)]; uw = Viw[3]
bra0w = prod_bra([vw, uw, vw, uw])
b0y = bra0w.dot(Ykn)
print('== <0_w|y> tensor ranks (site ranks all 1 <=> product over sites):', tensor_ranks(list(b0y / b0y[0])))
cw = None
idx = {w: i for i, w in enumerate(WORDS)}
cwv = [Fr(0)] * 256
for (y, x), val in M.items():
    cwv[idx[x]] += b0y[idx[y]] * val
print('expansion of <0_w|: nnz', sum(1 for x in cwv if x != 0), ' reproduces the bra:', list(np.array(cwv, dtype=object).dot(np.array(W, dtype=object))) == list(bra0w))
# eigenvalues of t^{Gc(w)} at the ket points, product formula
chw = Chain(G=Gc(Wt))
def ev(op, bra):
    r = bra.dot(op); p = next(i for i, x in enumerate(bra) if x != 0); e = r[p] / bra[p]
    return e, all(r[i] == e * bra[i] for i in range(256))
h, th = H0, TH0
f = {}
for s in range(4):
    for n in range(1, 4):
        if s in (0, 2) and n > 1: continue
        u = th[s] if s in (0, 2) else th[s] + n * h
        f[(s, n)] = ev(chw.t(1 if s in (0, 2) else n, u), bra0w)
print('eigenbra property at ket points:', all(v[1] for v in f.values()))
err = 0
for i, y in enumerate(WORDS):
    pr = Fr(1)
    for s in range(4):
        if y[s]: pr *= (f[(s, 1)][0] ** y[s]) if s in (0, 2) else f[(s, y[s])][0]
    err = max(err, abs(b0y[i] / (b0y[0] * pr) - 1))
print('<0_w|y> = <0_w|0bar> prod f: max rel err', err)
# z-independence of bra/ket words
Ww = sov_words(chw); print('bra words z-independent:', all(Ww[i][j] == W[i][j] for i in range(256) for j in range(256)))
