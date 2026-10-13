import numpy as np, pickle, sympy as sp, time
from sympy import prevprime
from core import *
from exact import crt_pair, ratrec_big
t0 = time.time()
z = z0
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
r = [l / lam[0] for l in lam]
b, bp = sp.symbols('b bp')
S = lambda x: sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else sp.Rational(x)
lamS = [S(x) for x in lam]; rS = [S(x) for x in r]
np_ = sp.Matrix([[1, 1, 1, rS[3]], [0, 1, 0, -rS[1]], [0, 0, 1, -rS[2]], [0, 0, 0, 1]])
Dp = sp.diag(1, rS[1], rS[2], 1)
g = np_ * Dp; gmT = g.inv().T
e = lambda i: sp.Matrix([[1 if k == i else 0 for k in range(1, 5)]])
def kron(vs):
    out = vs[0]
    for v in vs[1:]: out = sp.Matrix([[a * c for a in out for c in v]])
    return out
Wp = lambda a: [e(1) + a * e(4), e(4) - a * e(1)]
one = sp.Matrix([[1, 1, 1, 1]]); l = sp.Matrix([lamS])
lhs = kron([Wp(b)[0] * g, Wp(b)[1] * gmT, Wp(bp)[0] * g, Wp(bp)[1] * gmT])
Bdiag = kron([l + b * lamS[0] * e(4), one - b * e(1), l + bp * lamS[0] * e(4), one - bp * e(1)])
print('decomposition (W_b x W_bp) rho(n\'D\') == Bdiag/lam1^2 :', sp.expand(lhs - Bdiag / lamS[0] ** 2) == sp.zeros(1, 256))
V = sp.Matrix([[S(x) for x in row] for row in Vand(z)]); Vi = V.inv()
vz4 = V[:, 3].T; l1 = Vi[0, :]
Bcomp = kron([e(1) + b * lamS[0] * vz4, e(4) - b * l1, e(1) + bp * lamS[0] * vz4, e(4) - bp * l1])
Rinv = lambda vec: vec  # placeholder
# frame check: Bcomp R == Bdiag with R = V^-T x V x V^-T x V : check per pair factor: A: (e1 + b lam1 vz4) V^-T == l + b lam1 e4 ; F: (e4 - b l1) V == one - b e1
print('frame map A:', sp.simplify((e(1) + b * lamS[0] * vz4) * Vi.T - (l + b * lamS[0] * e(4))) == sp.zeros(1, 4),
      ' F:', sp.simplify((e(4) - b * l1) * V - (one - b * e(1))) == sp.zeros(1, 4))
print("n' =", np_.tolist(), " D' = diag", [1, rS[1], rS[2], 1])
# exact SoV coefficients of B(b): pieces via multi-prime
primes = []; q = 2 ** 26
while len(primes) < 20:
    q = prevprime(q); primes.append(q)
try:
    Ys = pickle.load(open('Y_mod20.pkl', 'rb'))
except Exception:
    Ys = {}
    for pp in primes:
        Rm = Ring('mod', pp); M, _ = sov_words(Chain(Rm)); Ys[pp] = Rm.inv(M)
    pickle.dump(Ys, open('Y_mod20.pkl', 'wb'))
print('Y mod 20 primes ready', time.time() - t0)
def piece_vec(bval_terms):
    pass
# pieces: pair1 vectors as polynomials in b: A = e1 + b A1, F = e4 + b F1 ; B = (e1 x e4) + b (A1 x e4 + e1 x F1) + b^2 (A1 x F1), all x (e1 x e4)
A1 = [lam[0] * x for x in Vand(z)[0:4] and [Fr(1)]]  # dummy (replaced below)
vz4F = [z[3] ** 3, z[3] ** 2, z[3], Fr(1)]
Vi_fr = [[Fr(int(sp.fraction(Vi[i, j])[0]), int(sp.fraction(Vi[i, j])[1])) for j in range(4)] for i in range(4)]
A1 = [lam[0] * x for x in vz4F]; F1 = [-x for x in Vi_fr[0]]
E1 = [1, 0, 0, 0]; E4 = [0, 0, 0, 1]
pieces_sites = {0: [[E1, E4]], 1: [[A1, E4], [E1, F1]], 2: [[A1, F1]]}
res = {}
for k, terms in pieces_sites.items():
    per = {}
    for pp in primes:
        Rm = Ring('mod', pp); acc = np.zeros(256, dtype=np.int64)
        for (a, f) in terms:
            bvec = kronvecs(Rm, [Rm.arr(a), Rm.arr(f), Rm.arr(E1), Rm.arr(E4)])
            acc = (acc + Rm.mm(bvec[None, :], Ys[pp])[0]) % pp
        per[pp] = acc
    vals = {}
    for x in range(256):
        if all(per[pp][x] == 0 for pp in primes[:3]): continue
        a_, m_ = int(per[primes[0]][x]), primes[0]
        for pp in primes[1:-1]: a_, m_ = crt_pair(a_, m_, int(per[pp][x]), pp)
        v = ratrec_big(a_, m_)
        ok = v is not None and (v.numerator - v.denominator * int(per[primes[-1]][x])) % primes[-1] == 0
        vals[wordsIdx[x]] = v if ok else None
    res[k] = vals
    print('b^%d piece: %d words, reconstruction failures %d' % (k, len(vals), sum(v is None for v in vals.values())))
supp = sorted(set().union(*[set(v.keys()) for v in res.values()]))
print('support size', len(supp), ' a3 values', sorted(set(x[2] for x in supp)), ' k4 values', sorted(set(x[3] for x in supp)))
pickle.dump(res, open('phys_pieces_exact.pkl', 'wb'))
print('| label (a₁,k₂,a₃,k₄) | coefficient of b¹ | coefficient of b² |'); print('|---|---|---|')
for x in supp:
    if x == (0, 0, 0, 0) and False: continue
    c1v = res[1].get(x, 0); c2v = res[2].get(x, 0)
    f = lambda v: '0' if v == 0 else ('?' if v is None else '%.6g' % float(v))
    print(f'| {x} | {f(c1v)} | {f(c2v)} |')
print('time', time.time() - t0)
