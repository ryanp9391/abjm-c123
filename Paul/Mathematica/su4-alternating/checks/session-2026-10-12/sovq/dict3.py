"""Full verification of the SoV <-> Q dictionary over all 256 words x 256 states.
Q convention ('mine'): sum_a (-1)^a tau_a(u) Q(u - a h) = 0, Q_i(u) = z_i^{u/h} Phi4(u) q_i(u), Phi4(u) = V4(u) Phi4(u-h).
RV/coordinator convention: Q^RV(u) = Q^mine(u - h).  Regularised values at site s: drop the word-independent factor
z_i^{th_s/h} Phi4(th_s - J h):  Qt_i(th_s - j h) = z_i^{-j} prod_{b=j}^{J-1} V4(th_s - b h) q_i(th_s - j h)."""
import pickle, itertools
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
from sov_l4 import WORDS
J = 5
def Qt(n, i, ths, j):
    return z[i - 1] ** (-j) * mp.fprod([V(4, ths - b * h) for b in range(j, J)]) * peval(allQ[n][(i,)], ths - j * h)
# --- direct statements ---
# (1) F site, any j:  tau1(thF) = z_j V4(thF) q_j(thF)/q_j(thF-h) = Q_j(thF)/Q_j(thF-h)   (N = 1)
mx = 0
for s in (1, 3):
    for n in range(256):
        for j in range(1, 5):
            pred = Qt(n, j, th[s], 0) / Qt(n, j, th[s], 1)
            mx = max(mx, abs(states[n].tau(1, th[s]) / pred - 1))
print('F site tau1(thF) = Q_j(thF)/Q_j(thF-h), all j, all states: max |ratio-1| =', mp.nstr(mx, 3))
# (2) A site, any triple I: tau_a(thA-h) = D_a/D_0, D_a = det_{i in I, j in {1,2,3,4}\{a+1}} Q_i(thA - j h)
def Dg(n, I, ths, pts):
    return mp.det(mp.matrix([[Qt(n, i, ths, j) for j in pts] for i in I]))
mx = 0
for s in (0, 2):
    for n in range(256):
        for I in itertools.combinations((1, 2, 3, 4), 3):
            D0 = Dg(n, I, th[s], (2, 3, 4))
            for a in (1, 2, 3):
                pts = [j for j in (1, 2, 3, 4) if j != a + 1]
                mx = max(mx, abs(states[n].tau(a, th[s] - h) * D0 / Dg(n, I, th[s], pts) - 1))
print('A site tau_a(thA-h) = D_a/D_0, all triples, all states: max |ratio-1| =', mp.nstr(mx, 3))
# (3) RV full product formula over GT data, for several sigma (Q_sigma(1), Q_sigma(2), Q_sigma(3)):
#     <Lambda|Psi>/<0|Psi> = prod_s prod_{k=1}^{3} det_{i,j<=k} Q^RV_{sigma(i)}(x^s_{kj}) / (same at Lambda = 0)
#     F site (RV x): x_k1 in {thF, thF+h} (k of the three dynamical ones at thF+h for word k: x_11 >= x_21 >= x_31),
#       x22 = x32 = thF - h, x33 = thF - 2h.   A site: x11 = x21 = thA, x22 = thA - h, x3j = (thA..thA-3h) minus thA - a h.
#     In 'mine' Q-argument: subtract h; in units of j = (th_s - argument)/h.
def xF(k):  # returns list over k_level of lists of j-values (mine: arg = thF - j h)
    mu = [1] * k + [0] * (3 - k)           # mu_11 >= mu_21 >= mu_31
    x1 = [mu[0]]; x2 = [mu[1], -1]; x3 = [mu[2], -1, -2]   # RV: thF + h*x
    return [[1 - v for v in x1], [1 - v for v in x2], [1 - v for v in x3]]   # mine: arg = thF + h*x - h
def xA(a):
    x1 = [0]; x2 = [0, -1]; x3 = [v for v in (0, -1, -2, -3) if v != -a]
    return [[1 - v for v in x1], [1 - v for v in x2], [1 - v for v in x3]]
def psi_site(n, sigma, ths, xs):
    p = mp.mpf(1)
    for k, js in enumerate(xs, start=1):
        p *= mp.det(mp.matrix([[Qt(n, sigma[i], ths, j) for j in js] for i in range(k)]))
    return p
mx = 0
for sigma in [(1, 2, 3), (4, 3, 2), (2, 4, 1), (3, 1, 4)]:
    for n in range(256):
        base = {}
        vals = {}
        for s in range(4):
            for lab in range(4):
                xs = xA(lab) if s in (0, 2) else xF(lab)
                vals[(s, lab)] = psi_site(n, sigma, th[s], xs)
        for w in WORDS:
            pred = mp.fprod([vals[(s, w[s])] / vals[(s, 0)] for s in range(4)])
            direct = states[n].tau(w[0], th[0] - h) * states[n].tau(1, th[1]) ** w[1] * states[n].tau(w[2], th[2] - h) * states[n].tau(1, th[3]) ** w[3]
            mx = max(mx, abs(direct / pred - 1))
    print('RV product formula, sigma =', sigma, ': max |ratio-1| over 256 words x 256 states =', mp.nstr(mx, 3))
