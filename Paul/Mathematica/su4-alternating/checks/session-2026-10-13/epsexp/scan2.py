"""scan2.py L [states] -- extended twisted scan: per (state, eps) store tau's, overlaps with W_a^L, <vac|J^k|Psi_eps>, T14 matrix elements,
odd-bridge chain data (k = 1, 2).  Usage: python scan2.py 2 ; python scan2.py 3 4,5,6,16,17,18,19,20,21,22 ; output scan2_L{L}.pkl"""
import numpy as np, sys, time, pickle, math
from ue2 import *
Lp = int(sys.argv[1]); states = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
phi = [0.9, 0.35, 0.35, -1.1]; thdir = [0.0, 0.37, -0.52]; c = 0.3
S = Setup(Lp)
a = (phi[0] - c) / (phi[0] - phi[3])
EPS = [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4]
if len(sys.argv) > 3: EPS = [float(x) for x in sys.argv[3].split(',')]
vac = S.vac()
VJ = [vac]
for k in range(2 * Lp + 1): VJ.append(S.J14.T @ VJ[-1])
rows = {}
t0 = time.time()
Us = [0.37, -0.29, 0.83]
for eps in EPS:
    ch, z = S.chain(eps, phi, thdir); y = z[1]
    kap = 1 - c * eps
    ap = a * (1 - kap * z[3]) / (1 - kap * z[0]); Nk = (1 - kap * z[0]) ** 2
    st = S.twisted_states(ch, states)
    braW = S.bra([a] * Lp); braL = S.bra([ap] + [a] * (Lp - 1))
    th = ch.th
    for i in sorted(st):
        psi, ov, sec, gap = st[i]
        taus, eres = S.pair_taus(ch, psi, sec)
        K = S.Kval(ch, taus, kap, y)
        w = braW @ psi; lhs = braL @ psi
        r = dict(eps=eps, ov=ov, gap=gap, eres=eres, K=K, w=w, lhs=lhs, N=Nk, taus=taus, rhs=K * w / Nk, psi=psi.astype(np.complex64) if False else None)
        r['vacJ'] = [VJ[k] @ psi for k in range(2 * Lp + 2)]
        r['tau1th'] = [tau_ket(ch, psi, 1, th[k]) for k in range(2 * Lp)]
        # exact identity (tau_vac(u) - tau_Psi(u)) <W_a|Psi> = a (z1-z4) <W_a|T14(u)|Psi>
        idn = []
        for u in Us:
            tv = (ch.apply_t(1, u, vac.reshape((4,) * ch.N)).reshape(-1) @ vac) / (vac @ vac)   # vacuum eigenvalue (diagonal twist): <vac|t = tv <vac|
            tp = tau_ket(ch, psi, 1, u)
            T14 = apply_t1G(ch, u, braW.reshape((4,) * ch.N), e41()).reshape(-1) @ psi
            T14v = apply_t1G(ch, u, vac.reshape((4,) * ch.N), e41()).reshape(-1) @ psi   # <vac|T14(u)|Psi>
            idn.append((u, tv, tp, T14, T14v))
        r['ident'] = idn
        # odd bridge chain, k = 1, 2
        odd = {}
        for k in (1, 2):
            bra = chain_bra(Lp, k, a, z, eps, c)
            nrm = eps * (c + (phi[3] if k % 2 == 0 else phi[0]))
            rhs_tau = chain_word_eig(Lp, k, th, r['tau1th'], c, eps) * w
            odd[k] = dict(lhs=bra @ psi, nrm=nrm, rhs=rhs_tau, wordeig=chain_word_eig(Lp, k, th, r['tau1th'], c, eps))
        r['odd'] = odd
        rows.setdefault(i, []).append(r)
    print('eps %.0e done (%.0f s)' % (eps, time.time() - t0)); sys.stdout.flush()
pickle.dump(dict(rows=rows, phi=phi, c=c, a=a, thdir=thdir, EPS=EPS), open('scan2_L%d.pkl' % Lp, 'wb'))
