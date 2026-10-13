# Task 1h: right-end pair words (pair L): Y1 = t1(thA_L), Y2 = t2(thA_L), Y3 = t3^(thA_L + h).
# (i) which products localise on pair L; (ii) intersection of their span with W_a^{L-1} (x) {W_b}; (iii) left x right products:
# intersection of span{w^(1) v^(L)} with {W_b (x) W_a^{L-2} (x) W_c}.
import numpy as np, itertools, sys
from fractions import Fraction as Fr
import eng, lp, modlin
from eng import P
from supp import support_pairs
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
Lp = int(sys.argv[1]) if len(sys.argv) > 1 else 3
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9), Fr(5, 6)][:Lp]
ch = lp.chain_abjm(Lp, thF, h, z, 'p')
a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
tA = ch.th[2 * Lp - 2]
Y = {'Y1': lambda v: ch.apply_t(1, tA, v), 'Y2': lambda v: ch.apply_t(2, tA, v), 'Y3': lambda v: ch.apply_t(3, tA + h, v, reg=True)}
words = {'1': ref}
for n in Y: words[n] = Y[n](ref)
for n1, n2 in itertools.combinations_with_replacement(['Y1', 'Y2', 'Y3'], 2):
    words[n1 + n2] = Y[n2](words[n1])
loc = {}
for n, v in words.items():
    S = support_pairs(ch, v, [a] * Lp)
    print('  %-6s changed pairs %s' % (n, S))
    if S is not None and set(S) <= {Lp - 1}: loc[n] = v
e = np.eye(4, dtype=np.int64)
u = [np.multiply.outer(e[0], e[3]), (np.multiply.outer(e[3], e[3]) - np.multiply.outer(e[0], e[0])) % P, (-np.multiply.outer(e[3], e[0])) % P]
wa = eng.prod_state(eng.W(a), 'p')
restL = wa
for _ in range(Lp - 2): restL = np.multiply.outer(restL, wa) % P
Ulast = [np.multiply.outer(np.multiply.outer(restL, wa) % P if False else restL, u[i]).reshape(-1) % P for i in range(3)] if False else None
# W_a^{L-1} (x) u_i
WL1 = wa
for _ in range(Lp - 2): WL1 = np.multiply.outer(WL1, wa) % P
U = [np.multiply.outer(WL1, u[i]).reshape(-1) % P for i in range(3)]
cols = [v.reshape(-1) for v in loc.values()] + [(-x) % P for x in U]
rk, cons, sol, null = modlin.solve(cols, np.zeros_like(cols[0]))
print('right-end localised words:', list(loc), ' rank', len(loc) - 0, '; intersection with W_a^{L-1} (x) span{u0,u1,u2}: dim', len(null))
# left x right
left = {'1': ref}
for i, x in enumerate(lp.Xs(ch, ref, 1)): left[['X1', 'X2', 'X3', 'X3X1'][i]] = x
mons = []
for nl, vl in left.items():
    for nr in loc:
        v = vl
        if nr != '1':
            for part in [nr[i:i + 2] for i in range(0, len(nr), 2)]: v = Y[part](v)
        mons.append(v.reshape(-1))
mid = np.ones(1, dtype=np.int64)
for _ in range(Lp - 2): mid = np.multiply.outer(mid, wa) % P
mid = mid.reshape(-1)
U2 = [np.multiply.outer(np.multiply.outer(u[i].reshape(-1), mid) % P, u[j].reshape(-1)).reshape(-1) % P for i in range(3) for j in range(3)]
cols = mons + [(-x) % P for x in U2]
rk, cons, sol, null = modlin.solve(cols, np.zeros_like(cols[0]))
print('left x right: %d monomials, rank of [monomials | 9 tensors u_i (x) W_a^{L-2} (x) u_j] = %d, intersection dim = %d' % (len(mons), rk, len(null)))
for nv in null:
    print('    d_ij:', [str(modlin.ratrec(int(x))) for x in nv[len(mons):]])
