import time, itertools
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
from lib import *
from words import factor_list
mp.mp.dps = 30
L = 4; a0 = Fr(1, 3)
chx = build_chain(L, dt=Fraction)      # exact rational arithmetic (stronger than 30 digits)
th = TH6; h = H; z = Z4
def kron_obj(vs):
    out = np.array([Fr(1)], dtype=object)
    for v in vs: out = np.array([x * y for x in out for y in v], dtype=object)
    return out
def Wp(a_list):
    vs = []
    for a in a_list:
        vs += [[Fr(1), 0, 0, a], [-a, 0, 0, Fr(1)]]
    return kron_obj(vs)
psi0 = Wp([a0, a0])
F = factor_list(chx)
t0 = time.time()
single = [chx.apply(f['a'], f['u'], psi0) for f in F]
print('exact singles built', time.time() - t0, flush=True)
# targets (exact)
T = {}
for a1 in (Fr(-2, 3), Fr(0), Fr(2, 3), Fr(5, 4)): T[f"T1 a'={a1}"] = Wp([a1, a0])
v = [Fr(1), 0, 0, a0]; u = [-a0, 0, 0, Fr(1)]; rest = kron_obj([v, u])
for k in (Fr(1, 3), Fr(-2, 5), Fr(5, 3)):
    vD = [v[i] * (1 - k * z[i]) for i in range(4)]; uD = [u[i] / (1 - k * z[i]) for i in range(4)]
    T[f"T2 k={k}"] = kron_obj([vD, uD]) if False else np.array([x * y for x in kron_obj([vD, uD]) for y in rest], dtype=object)
for k in (Fr(1, 3),):
    vD = [v[i] * (1 - k / z[i]) for i in range(4)]
    T[f"T3a k={k}"] = np.array([x * y for x in kron_obj([vD, u]) for y in rest], dtype=object)

from modlib import ModChain, PRIMES
from modrun import Wex_pair, vec_mod
from pipeline import build_positive_words
from cyc import ModBasis
p = PRIMES[0]
mc = ModChain(L, TH6, H, Z4, p)
v0, u0 = Wex_pair(Fraction(1, 3))
psim = vec_mod([v0, u0] * 2, p)
Fm = factor_list(mc)
ws, _, _ = build_positive_words(mc, psim, Fm, ('W1', 'W2'), verbose=False)
Wm = ws.matrix(); labs = ws.labels
def independent_words(maxf):
    B = ModBasis(mc.dim, p); sel = []
    for j in sorted(range(len(labs)), key=lambda j: len(labs[j])):
        if len(labs[j]) > maxf: continue
        if B.add(Wm[:, [j]].T.copy())[0]: sel.append(j)
    return sel
def exact_word(lab):
    vec = psi0
    for (i, e) in lab: vec = chx.apply(F[i]['a'], F[i]['u'], vec)
    return vec
def to_mp(vec): return [mp.mpf(x.numerator) / mp.mpf(x.denominator) if isinstance(x, Fr) else mp.mpf(x) for x in vec]
def dot(a, b): return mp.fsum(x * y for x, y in zip(a, b))
def resid(cols, b):
    """modified Gram-Schmidt (twice) at 30 digits; residual of b after projecting on span(cols), relative"""
    Q = []
    for c in cols:
        w = to_mp(c)
        for _ in range(2):
            for q in Q:
                cf = dot(q, w); w = [wi - cf * qi for wi, qi in zip(w, q)]
        nw = mp.sqrt(dot(w, w))
        if nw < mp.mpf(10) ** (-20) * mp.sqrt(dot(to_mp(c), to_mp(c))): continue   # dependent column
        Q.append([wi / nw for wi in w])
    w = to_mp(b); nb = mp.sqrt(dot(w, w))
    for _ in range(2):
        for q in Q:
            cf = dot(q, w); w = [wi - cf * qi for wi, qi in zip(w, q)]
    return mp.sqrt(dot(w, w)) / nb, len(Q)
for maxf in (1, 2):
    sel = independent_words(maxf)
    cols = [exact_word(labs[j]) for j in sel]
    print(f'<= {maxf} factors: independent words (mod p pivot selection) {len(sel)}; 30-digit QR least squares on the exact bras', flush=True)
    for k, tv in T.items():
        rr, nq = resid(cols, tv); print('   ', k, 'relative residual', mp.nstr(rr, 5), '(independent columns used: %d)' % nq, flush=True)
