import sys
from fractions import Fraction as Fr
import numpy as np
from modlib import *
from modrun import mod_targets, Wex_pair, vec_mod
from cyc import ModBasis
from words import factor_list
def run(L, th, p=PRIMES[0]):
    mc = ModChain(L, th, H, Z4, p)
    a0 = Fr(1, 3); v, u = Wex_pair(a0)
    psi0 = vec_mod([v, u] * (L // 2), p)
    T, w0 = mod_targets(L, a0, Z4, p)
    t3 = lambda x: mc.apply(3, th[0] - H, x)
    t1b = lambda x: mc.apply(1, th[1], x)
    w_1 = psi0; w_t3 = t3(psi0); w_t3t1 = t1b(w_t3); w_t1 = t1b(psi0)
    rest = [Wex_pair(a0)[i] for _ in range(L // 2 - 1) for i in (0, 1)]
    Sfull = np.array([sum(0 for _ in [])], dtype=float)
    Sv = []
    out = [Fr(0)] * 16
    for i in range(4): out[i * 4 + i] = Fr(1)
    kr = out
    for vs in rest: kr = [x * y for x in kr for y in vs]
    S = np.array([int(Fr(x).numerator % p) * pow(int(Fr(x).denominator % p), -1, p) % p for x in kr], dtype=float)
    B4 = ModBasis(mc.dim, p); B4.add(np.stack([w_1, w_t3, w_t3t1, w_t1]))
    B5 = ModBasis(mc.dim, p); B5.add(np.stack([w_1, w_t3, w_t3t1, w_t1, S]))
    r4 = {k: B4.contains(v_) for k, v_ in T.items() if k.startswith('T1') or k.startswith('T2')}
    r5 = {k: B5.contains(v_) for k, v_ in T.items() if k.startswith('T1') or k.startswith('T2')}
    print(f'L sites {L}: rank(4 words) {B4.rank}, rank(4 words + S) {B5.rank}; S in span(4 words): {B4.contains(S)}')
    print('   T1/T2 in span of 4 localised words:', all(r4.values()), ' in span of 4 words + S(x)w0:', all(r5.values()))
for L in (4, 6):
    for th in (TH6, TH6b, TH6c):
        run(L, th)
