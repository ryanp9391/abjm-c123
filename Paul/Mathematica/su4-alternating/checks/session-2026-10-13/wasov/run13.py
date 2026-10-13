"""Even bridge M = 1 at L = 2 with <W_a|Psi> replaced by the Q-closed form (F-only SoV bra, bidegree (2,2) in rho2, rho4):
   w_Q = <0|Psi> (a/kappa0)^n (ut.P(rho2,rho4).ut)/norm  [diag frame: <0| = (lam x 1)^2],  C_Q(eps) = K_kappa(tau) w_Q/(1-kappa z1)^2.
z2 != z3 (phi2 != phi3) is needed for <0|Psi> (companion frame); K_kappa uses y = (z2+z3)/2 (lpairs Task 3c)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'epsexp')))
from ue_base import *
Lp = 2
S = Setup(Lp)
phi = [0.9, 0.35, 0.6, -1.1]; thdir = [0.0, 0.37]; c = 0.3
EPSL = [float(x) for x in os.environ.get("EPSL", "1e-2,5e-3,2e-3").split(",")]
if len(sys.argv) > 1: phi = [float(x) for x in sys.argv[1].split(',')]
a = (phi[0] - c) / (phi[0] - phi[3])
def Pab(nu, hh, r2, r4):
    return np.array([[(nu / (nu + hh)) ** 2 * r2 ** 2 * r4 ** 2, nu / (nu + hh) * r2 ** 2 * r4, (nu - hh) / nu * r2 ** 2 + hh / (nu + hh) * r2 * r4],
                     [nu / (nu + hh) * r2 * r4 ** 2 - hh / (nu + hh) * r2 ** 2 * r4, nu / (nu + hh) * r2 * r4, r2],
                     [r4 ** 2 - hh / nu * (1 + hh / (nu + hh)) * r2 * r4 - hh / nu * (nu - hh) / nu * r2 ** 2, r4 - hh / nu * r2, 1.0]], dtype=complex)
print('phi', phi, 'c', c, 'a', a)
res = {}
for eps in EPSL:
    ch, z = S.chain(eps, phi, thdir)
    y = 0.5 * (z[1] + z[2])
    kap = 1 - c * eps; Nk = (1 - kap * z[0]) ** 2
    th = ch.th
    nu = th[1] - th[3]
    lam = np.array([1 / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)])
    bra0 = eng.prod_state([lam, np.ones(4), lam, np.ones(4)], 'c').reshape(-1)
    kap0 = -((z[0] - z[1]) * (z[0] - z[2])) / ((z[3] - z[1]) * (z[3] - z[2]))
    norm = lam[0] ** 2 * ((z[3] - z[1]) * (z[3] - z[2])) ** 2
    ut = np.array([1.0, -(z[1] + z[2]), z[1] * z[2]])
    st = S.twisted_states(ch)
    braW = S.bra([a] * Lp)
    for i in sorted(st):
        psi, ov, sec, gap = st[i]
        n = nval(Lp, S.DYN[i])
        if n != 2: continue
        taus, eres = S.pair_taus(ch, psi, sec)
        K = S.Kval(ch, taus, kap, y)
        # t1(theta2), t1(theta4) eigenvalues
        r2, _ = S.eigval(ch, psi, sec, 1, th[1]); r4, _ = S.eigval(ch, psi, sec, 1, th[3])
        rho2 = r2 / (h * (th[0] - th[1]) * (th[1] - th[2]) * (nu - h)); rho4 = r4 / (h * (th[0] - th[3]) * nu * (th[2] - th[3]))
        P = Pab(nu, h, rho2, rho4)
        wQ = (bra0 @ psi) * (a / kap0) ** n * (ut @ P @ ut) / norm
        w = braW @ psi
        res.setdefault(i, []).append((eps, K * w / Nk, K * wQ / Nk, abs(w - wQ) / abs(w), K, w, wQ))
        print('eps %.0e state #%d %s: |w-wQ|/|w| = %.2e   C_direct = %+.6f%+.6fi   C_Q = %+.6f%+.6fi   (data %+.6f)' % (
            eps, i, S.DYN[i], abs(w - wQ) / abs(w), (K * w / Nk).real, (K * w / Nk).imag, (K * wQ / Nk).real, (K * wQ / Nk).imag, S.C0[i][2].real), flush=True)
print('\nRichardson (eps 1e-2,5e-3,2e-3, linear + quadratic fit) of C_Q and C_direct, data = C0[i][2]:')
for i, rr in res.items():
    e = np.array([r[0] for r in rr])
    for name, col in (('direct', 1), ('Q', 2)):
        v = np.array([r[col] for r in rr])
        co = np.polyfit(e, v, 2)
        print('  state #%d %s C_%s(0) = %+.6f%+.6fi   (data %+.6f)' % (i, S.DYN[i], name, co[-1].real, co[-1].imag, S.C0[i][2].real))
