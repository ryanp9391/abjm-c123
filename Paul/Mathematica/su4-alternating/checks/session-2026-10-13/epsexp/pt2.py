"""pt2.py L states -- exact second-order formula for n = 2 states in terms of T14 form factors:
  omega_2 (tau_vac(u) - tau_Psi0(u)) = 2 (phi1-phi4)^2 sum_m <vac|T14(u)|m> <m~|T14(u1)|Psi0> / (tau_m(u1) - tau_Psi0(u1)),
  m = eigenvectors of t1(u1) (untwisted) in the sector of J Psi0; and  <vac|T14(u)|m0> = kappa_m (tau_vac(u) - tau_m(u)) for n = 1 states."""
import numpy as np, sys, pickle
from ue2 import *
np.set_printoptions(linewidth=220, precision=6, suppress=True)
Lp = int(sys.argv[1]); states = [int(x) for x in sys.argv[2].split(',')]
phi = [0.9, 0.35, 0.35, -1.1]
S = Setup(Lp); N = S.N; D = S.D; vac = S.vac()
ch0, z0 = S.chain(0.0, phi, [0, .37, -.52])
u1 = 0.41; us = [0.37, -0.29, 0.83]
def T14mat(u, secA, secB, chunk=64):
    out = np.zeros((len(secA), len(secB)), dtype=complex)
    for b0 in range(0, len(secA), chunk):
        idx = secA[b0:b0 + chunk]; B = len(idx)
        E = np.zeros((D, B), dtype=complex); E[idx, np.arange(B)] = 1
        res = apply_t1G(ch0, u, E.reshape((4,) * N + (B,)), e41()).reshape(D, B)   # res[q,b] = <e_p|T14|e_q>
        out[b0:b0 + B] = res[secB].T
    return out
def tauvac(u): return (ch0.apply_t(1, u, vac.reshape((4,) * N)).reshape(-1) @ vac) / (vac @ vac)
out = {}
for i in states:
    n = nval(Lp, S.DYN[i]); w = tuple(S.W0[i])
    secB = S.SECT[w]; wA = (w[0] + 1, w[1], w[2], w[3] - 1); secA = S.SECT[wA]
    psi0 = S.PSI0[i][secB]
    tpsi = {u: tau_ket(ch0, S.PSI0[i], 1, u) for u in us + [u1]}
    M1 = S.sector_matrix(ch0, [(1.0, 1, u1, False)], secA)
    ev, R = np.linalg.eig(M1); Rinv = np.linalg.inv(R)
    T1 = T14mat(u1, secA, secB); vecA = T1 @ psi0
    comp = Rinv @ vecA
    bu = {u: apply_t1G(ch0, u, vac.reshape((4,) * N), e41()).reshape(-1)[secA] for u in us}
    # X_1(u) from the PT sum
    print('#%d %s n=%s  dim(sec)=%d dim(J sec)=%d' % (i, S.DYN[i], n, len(secB), len(secA)))
    for u in us:
        Dt = tauvac(u) - tpsi[u]
        terms = []
        for m in range(len(ev)):
            c_m = comp[m]
            if abs(c_m) < 1e-9: continue
            term = (bu[u] @ R[:, m]) * c_m / (ev[m] - tpsi[u1])
            if abs(term) > 1e-9: terms.append((ev[m], term, ((bu[u] @ R[:, m]))))
        rhs = 2 * (phi[0] - phi[3]) ** 2 * sum(t[1] for t in terms) / Dt
        print('   u=%.2f  omega_2 predicted = %s   (%d contributing eigenvectors; tau_m(u1)=%s)' % (u, np.round(rhs, 6), len(terms), [np.round(t[0], 4) for t in terms]))
    # C_{M=1} prediction: C = K1 * omega_2 / (2 (phi1-phi4)^2) = K1/Dtau(u) * sum_m <vac|T14(u)|m><m~|T14(u1)|Psi0>/(tau_m(u1)-tau_Psi0(u1))
    K1 = pickle.load(open('scan_L%d_c0.3.pkl' % Lp, 'rb'))['K1'][i][0]
    u = us[0]; Dt = tauvac(u) - tpsi[u]
    ssum = sum((bu[u] @ R[:, m]) * comp[m] / (ev[m] - tpsi[u1]) for m in range(len(ev)))
    Cpred = K1 * ssum / Dt
    # Casimir of contributing eigenvectors
    Cm = None
    for ii in range(4):
        for jj in range(4):
            t_ = S.global_gen(ii, jj) @ S.global_gen(jj, ii); Cm = t_ if Cm is None else Cm + t_
    Cs = Cm[secA][:, secA].toarray(); casm = np.real(np.diag(Rinv @ Cs @ R))
    contrib = [(round(float(casm[m]), 3)) for m in range(len(ev)) if abs((bu[u] @ R[:, m]) * comp[m]) > 1e-9]
    print('   C_{M=1} predicted = K1 * sum/Dtau = %s ; data C(M=1) = %s ; Casimirs of contributing eigenvectors: %s' % (np.round(Cpred, 6), np.round(S.C0[i][2 * Lp - 2], 6), sorted(set(contrib))))
    out[i] = dict(n=n)
