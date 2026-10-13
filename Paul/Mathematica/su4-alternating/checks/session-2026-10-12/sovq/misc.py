import pickle, itertools
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
from qfun import *
from poly import *
from sov_l4 import WORDS, TH0, H0, Z0
import frames
data = pickle.load(open('eigen.pkl', 'rb')); states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
# (c) polynomial-gauge normalisation: N_x = prod_A prod_{b=1}^{a} V4(thA - b h) * prod_F V4(thF)^k
#     <x|Psi>/<0|Psi> = N_x * prod_F [qh_j(thF)/qh_j(thF-h)]^k * prod_A det[qh_i(thA - j h)]_{j in S_a}/det[...]_{S_0}, qh_i(u) = z_i^{u/h} q_i(u)
def qh(n, i, u, ref): return z[i - 1] ** ((u - ref) / h) * peval(allQ[n][(i,)], u)
mx = 0
for n in range(256):
    for w in WORDS:
        val = mp.mpf(1)
        for s in (0, 2):
            a = w[s]
            Sa = [j for j in (1, 2, 3, 4) if j != a + 1]
            Na = mp.fprod([V(4, th[s] - b * h) for b in range(1, a + 1)])
            num = mp.det(mp.matrix([[qh(n, i, th[s] - j * h, th[s]) for j in Sa] for i in (1, 2, 3)]))
            den = mp.det(mp.matrix([[qh(n, i, th[s] - j * h, th[s]) for j in (2, 3, 4)] for i in (1, 2, 3)]))
            val *= Na * num / den
        for s in (1, 3):
            val *= (V(4, th[s]) * qh(n, 4, th[s], th[s]) / qh(n, 4, th[s] - h, th[s])) ** w[s]
        direct = states[n].tau(w[0], th[0] - h) * states[n].tau(1, th[1]) ** w[1] * states[n].tau(w[2], th[2] - h) * states[n].tau(1, th[3]) ** w[3]
        mx = max(mx, abs(val / direct - 1))
print('(c) polynomial-gauge dictionary with N_x: max |ratio-1| =', mp.nstr(mx, 3))
thF = [Fr(t) for t in TH0]; hF = H0
def V4F(u): return (u - thF[0]) * (u - thF[2]) * (u - thF[1] - hF) * (u - thF[3] - hF)
print('   N_F(k=1): V4(th2) =', V4F(thF[1]), ' V4(th4) =', V4F(thF[3]))
for s in (0, 2):
    print('   N_A site %d (a=1,2,3):' % (s + 1), [str(math.prod([V4F(thF[s] - b * hF) for b in range(1, a + 1)])) for a in (1, 2, 3)])
# (f) <0|Psi> for normalised eigenkets
ov = []
we = pickle.load(open('words_exact.pkl', 'rb')); W0 = we['words'][0]
for n in range(256):
    k = frames.companion_ket(n)
    nrm = mp.sqrt(mp.fsum(abs(c) ** 2 for c in k))
    ov.append(abs(mp.fsum(frames.tomp(W0[c]) * k[c] for c in range(256) if W0[c] != 0)) / nrm)
print('(f) |<0|Psi>|/|Psi| over 256 states: min', mp.nstr(min(ov), 5), 'max', mp.nstr(max(ov), 5), '# below 1e-30:', sum(1 for o in ov if o < 1e-30))
pickle.dump(ov, open('vacov.pkl', 'wb'))
# (d) diagonal of the Gram matrix: factorisation over sites?
ms = pickle.load(open('measure_exact.pkl', 'rb')); Gb = ms['Gb']
D = {x: Gb[(x, x)] for x in WORDS}
# additivity test of log|D| : D(x) D(x') = D(x'') D(x''') for site-label swaps
bad = 0; tot = 0
for x in WORDS:
    for xp in WORDS:
        for s in range(4):
            u1 = list(x); u2 = list(xp); u1[s], u2[s] = xp[s], x[s]
            u1, u2 = tuple(u1), tuple(u2)
            tot += 1
            if D[x] * D[xp] != D[u1] * D[u2]: bad += 1
print('(d) diagonal Gram entries factorise over sites (exact swap test):', bad == 0, '(%d violations of %d)' % (bad, tot))
