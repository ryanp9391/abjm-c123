"""untw.py -- untwisting limit of the rotated pair 1 (z = 1 + eps phi, w = 1 + eps omega, kappa = 1 - c eps).  Python fallback for notebook section 8."""
from helpers import *
import math
out = open('untw_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def nrm(x):
    n = math.sqrt(float(sum(a * a for a in x))); return [float(a) / n for a in x]
def dev(T, Tp):
    a = nrm(T); b = nrm(Tp)
    return min(math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))), math.sqrt(sum((x + y) ** 2 for x, y in zip(a, b))))
om = [Fr(7, 10), Fr(-3, 10), Fr(1, 5), Fr(1, 2)]; cc = Fr(4, 5)
P('generic phi = (9/10, 7/20, 7/10, -11/10) (phi2 != phi3): deviation of T_kappa = [v(1-kappa Lambda)] (x) [u (1-kappa Lambda)^-1] (v,u = diagonal-frame twisted vacuum at generic z) from (vhat D) (x) (uhat D^-1)')
phi = [Fr(9, 10), Fr(7, 20), Fr(7, 10), Fr(-11, 10)]
P('eps | D = c - phi | D = c + phi')
for e in (2, 3, 4, 5):
    ep = Fr(1, 10 ** e); zz = [1 + ep * p for p in phi]; ww = [1 + ep * o for o in om]; kap = 1 - cc * ep
    v, u = tw_form(ww, zz); T = kron([v[j] * (1 - kap * zz[j]) for j in range(4)], [u[j] / (1 - kap * zz[j]) for j in range(4)])
    vh = [prod_([(om[0] - phi[m]) / (phi[j] - phi[m]) for m in range(4) if m != j]) for j in range(4)]
    uh = [prod_([phi[j] - om[m] for m in range(3)]) for j in range(4)]
    devs = []
    for s in (-1, 1):
        Tp = kron([vh[j] * (cc + s * phi[j]) for j in range(4)], [uh[j] / (cc + s * phi[j]) for j in range(4)])
        devs.append(dev(T, Tp))
    P(e, '%.3e' % devs[0], '%.3e' % devs[1])
P('phi2 = phi3 = 7/20 (ABJM-point twist): reference = limit of the twisted vacuum (0,1,-1,0) (x) u_w; deviation of T_kappa from (0,1,-1,0) (x) (uhat_j/(c-phi_j))')
phi = [Fr(9, 10), Fr(7, 20), Fr(7, 20), Fr(-11, 10)]
for e in (2, 3, 4, 5):
    ep = Fr(1, 10 ** e); zz = [1 + ep * p for p in phi]; ww = [1 + ep * o for o in om]; kap = 1 - cc * ep
    u = tw_u(ww, zz); T = kron([Fr(0), 1 * (1 - kap * zz[1]), -(1 - kap * zz[2]), Fr(0)], [u[j] / (1 - kap * zz[j]) for j in range(4)])
    uh = [prod_([phi[j] - om[m] for m in range(3)]) for j in range(4)]
    Tp = kron([Fr(0), Fr(1), Fr(-1), Fr(0)], [uh[j] / (cc - phi[j]) for j in range(4)])
    P(e, '%.3e' % dev(T, Tp))
P('size of vhat components for phi2 -> phi3 (1/(phi2-phi3) pole): phi3 = 7/20 + d; vhat ~', )
for d in (Fr(1, 10), Fr(1, 100), Fr(1, 1000)):
    phi = [Fr(9, 10), Fr(7, 20), Fr(7, 20) + d, Fr(-11, 10)]
    vh = [prod_([(om[0] - phi[m]) / (phi[j] - phi[m]) for m in range(4) if m != j]) for j in range(4)]
    n = nrm(vh)
    P('d =', float(d), 'normalised vhat =', ['%.5f' % x for x in n])
out.close()
