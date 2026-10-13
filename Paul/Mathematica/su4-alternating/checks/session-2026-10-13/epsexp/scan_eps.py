"""scan_eps.py L -- orders in eps of w = <W_a^L|Psi_eps>, K = K_kappa(tau_Psi_eps), LHS = <W_a' W_a^{L-1}|Psi_eps>, N = (1-kappa z1)^2
for the linear prescription kappa = 1 - c eps, a = (phi1-c)/(phi1-phi4); plus K_1(tau_Psi0) at the exactly untwisted point,
and the limit formula C_{M=1} = K_1(tau_Psi0) * omega_n / (n! (phi1-phi4)^2) with omega_n = lim n! w/(a eps)^n (n = L-p)."""
import numpy as np, sys, time, pickle
from ue_base import *
Lp = int(sys.argv[1]) if len(sys.argv) > 1 else 2
phi = [0.9, 0.35, 0.35, -1.1]; thdir = [0.0, 0.37, -0.52]; c = 0.3
if len(sys.argv) > 2: phi = [float(x) for x in sys.argv[2].split(',')]
if len(sys.argv) > 3: c = float(sys.argv[3])
if len(sys.argv) > 4: thdir = [float(x) for x in sys.argv[4].split(',')]
S = Setup(Lp)
a = (phi[0] - c) / (phi[0] - phi[3])
EPS = [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4]
cidx = 2 * Lp - 2    # C index of the M = 1 even bridge: l = 2L-2 vacuum sites then one W_-1 pair = <vac^{L-1} W_-1| = <W_-1 vac^{L-1}| (cyclic)
rows = {}
t0 = time.time()
for eps in EPS:
    ch, z = S.chain(eps, phi, thdir); y = z[1]
    kap = 1 - c * eps
    ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
    Nk = (1 - kap * z[0]) ** 2
    st = S.twisted_states(ch)
    braW = S.bra([a] * Lp); braL = S.bra([ap] + [a] * (Lp - 1))
    for i in sorted(st):
        psi, ov, sec, gap = st[i]
        taus, eres = S.pair_taus(ch, psi, sec)
        K = S.Kval(ch, taus, kap, y)
        w = braW @ psi; lhs = braL @ psi
        rows.setdefault(i, []).append(dict(eps=eps, ov=ov, gap=gap, eres=eres, K=K, w=w, lhs=lhs, N=Nk, ap=ap, taus=taus,
                                          rhs=K * w / Nk))
    print('eps %.0e done (%.0f s)' % (eps, time.time() - t0)); sys.stdout.flush()
# untwisted K_1
ch0 = untwisted_chain(Lp)
K1 = {}
for i in range(len(S.PSI0)):
    sec = S.SECT[tuple(S.W0[i])]
    taus, eres = S.pair_taus(ch0, S.PSI0[i], sec)
    K1[i] = (K1_untwisted(ch0, taus), taus, eres)
pickle.dump(dict(rows=rows, K1=K1, phi=phi, c=c, a=a, thdir=thdir), open('scan_L%d_c%g.pkl' % (Lp, c), 'wb'))

def slope(xs, ys):
    xs = np.log(np.array(xs)); ys = np.log(np.abs(np.array(ys)) + 1e-300)
    return np.polyfit(xs, ys, 1)[0]

print('\n==== L = %d, phi = %s, c = %g, a = %.5f ====' % (Lp, phi, c, a))
print('state dyn n | slopes (log-log over eps 1e-2..1e-4): w, K, LHS | K_1(tau_Psi0) untwisted | lim K (eps=1e-4) | w/eps^n (1e-4) | C formula | C data | |RHS-LHS| max')
for i in sorted(rows):
    rr = rows[i]; n = nval(Lp, S.DYN[i])
    eps = [r['eps'] for r in rr]
    sw = slope(eps, [r['w'] for r in rr]); sK = slope(eps, [r['K'] for r in rr]); sL = slope(eps, [r['lhs'] for r in rr])
    k1 = K1[i][0]
    last = rr[-1]
    import math
    omega = last['w'] / (a * last['eps']) ** n * math.factorial(n) if n is not None else 0
    Cform = (k1 * omega / (math.factorial(n) * (phi[0] - phi[3]) ** n) if n == 2 else (0 if n != 0 else 1))
    dmax = max(abs(r['rhs'] - r['lhs']) for r in rr)
    print('#%2d %s n=%s | %+.2f %+.2f %+.2f | K1=%+.6f%+.6fi | K(1e-4)=%+.6f%+.6fi | omega_n=%+.6f%+.6fi | Cform=%+.6f%+.6fi | Cdata=%+.6f%+.6fi | %.0e'
          % (i, S.DYN[i], n, sw, sK, sL, k1.real, k1.imag, last['K'].real, last['K'].imag, omega.real, omega.imag,
             complex(Cform).real, complex(Cform).imag, S.C0[i][cidx].real, S.C0[i][cidx].imag, dmax))
print('\nK_1 residuals (eigen-residual at the untwisted point):', max(K1[i][2] for i in K1))
print('untwisted pair-point eigenvalues (tau_F, tau2(u0), tau3^(u0), x = tau1(u0)) and K_1:')
for i in range(len(S.PSI0)):
    t = K1[i][1]
    print('#%2d %s E=%7.3f  tauF=%+.6f%+.6fi  tau2=%+.6f%+.6fi  tau3h=%+.6f%+.6fi  x=%+.6f%+.6fi  K1=%+.6f%+.6fi'
          % (i, S.DYN[i], S.E0[i], t[0].real, t[0].imag, t[1].real, t[1].imag, t[2].real, t[2].imag, t[3].real, t[3].imag, K1[i][0].real, K1[i][0].imag))
