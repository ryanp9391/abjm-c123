# Task 4c: (1) factorisation of the per-pair eigenvalue polynomial K(kappa) = (1 - kappa rho1)(1 - kappa rho2), L = 2, z2 = z3;
#          (2) continuity: A-site Cramer formula tau_a(thA - h) at theta_A = theta_F + 2h + delta -> ABJM values.
import numpy as np, itertools
import eng, lp, qf, qdict
h = 0.5; thF = [-0.6, -1.375]
z = [1.5, 1.25, 1.25, 4 / 7]; y = z[1]
ch = lp.chain_abjm(2, thF, h, z, 'c')
V, ev, gaps = qf.joint_eigen(ch)
uF = thF[0]; u0 = uF + h
t1, _ = qf.eigvals_at(ch, V, 1, uF); t2, _ = qf.eigvals_at(ch, V, 2, u0); t3, _ = qf.eigvals_at(ch, V, 3, u0, reg=True)
l1, l2, l3 = eng.lam(ch, 1, uF, 1), eng.lam(ch, 2, u0, 1), eng.lam(ch, 3, u0, 1)
Sig = t1 / (h**2 * l1) + t2 / (2 * h**4 * y * l2) + t3 / (6 * h**5 * y**2 * l3)
Pi = t3 * t1 / (24 * h**7 * y**2 * l1 * l3)
# K = s^2 + k s Sig + k^2 Pi, s = 1 - k(z1+z4): K = 1 - k(2S - Sig) + k^2 (S^2 - S Sig + Pi), S = z1 + z4
S = z[0] + z[3]
b = 2 * S - Sig; c = S * S - S * Sig + Pi
disc = b * b - 4 * c
realstate = np.abs(t1.imag) + np.abs(t2.imag) + np.abs(t3.imag) < 1e-9
print('L=2, z2=z3: states with real (t1,t2,t3): %d of 256;  of these with negative discriminant (complex rho1,2): %d;  zero discriminant: %d'
      % (realstate.sum(), np.sum(realstate & (disc.real < -1e-9)), np.sum(np.abs(disc) < 1e-9)))
print('   vacuum-like examples (disc ~ 0) rho = b/2:', np.round(b[np.abs(disc) < 1e-9][:6] / 2, 6))
# distinct (Sigma, Pi) pairs:
print('   number of distinct K polynomials: %d' % len(set((round(x.real, 8), round(x.imag, 8), round(w.real, 8), round(w.imag, 8)) for x, w in zip(Sig, Pi))))
# (2) continuity of the A-site Cramer dictionary near the ABJM point (generic twist)
zg = [1.5, 1.25, 2 / 3, 4 / 7]
for delta in (1e-2, 1e-3, 1e-4, 0.0):
    th = [thF[0] + 2 * h + delta, thF[0], thF[1] + 2 * h, thF[1]]
    chd = eng.Chain(2, th, h, zg, 'c')
    Vd, evd, gd = qf.joint_eigen(chd)
    a1, _ = qf.eigvals_at(chd, Vd, 1, th[0] - h); a3, _ = qf.eigvals_at(chd, Vd, 3, th[0] - h)
    a3h, _ = qf.eigvals_at(chd, Vd, 3, th[0] - h, reg=True)
    # match states to delta = 0 later by sorting on the generic combination eigenvalue
    order = np.argsort(evd.real + 1e-3 * evd.imag)
    if delta == 1e-2: ref = None
    print('delta=%.0e: tau1(thA-h) first states %s ; tau3(thA-h)/delta %s ; tau3^ %s' % (delta, np.round(a1[order][:3], 6),
          np.round(a3[order][:3] / delta, 6) if delta else '-', np.round(a3h[order][:3], 6)))
