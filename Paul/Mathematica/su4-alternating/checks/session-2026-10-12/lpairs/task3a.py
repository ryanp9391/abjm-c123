# Task 3a: direction independence.  kappa = 1 - c eps exactly (linear), a fixed by a'-a = -1 at leading order:
#   a = (phi1 - c)/(phi1 - phi4)  [M = 1, target C4];  three directions phi (phi2 = phi3), two c.
import numpy as np, sys, time, pickle
import untw
thdir = [0.0, 0.37, -0.52]
PHIS = [(0.9, 0.35, 0.35, -1.1), (-0.7, 1.3, 1.3, 0.4), (1.5, -0.6, -0.6, 0.2)]
CS = (0.3, 2.0)
STATES = [0, 4, 16, 17, 18, 20, 22, 25]
out = []
for phi in PHIS:
    for eps in (1e-2, 1e-3, 1e-4):
        t0 = time.time()
        ch, z = untw.make_chain(eps, phi, thdir); y = z[1]
        st = untw.twisted_states(ch, STATES)
        for i in STATES:
            psi, ov, sec, gap = st[i]
            taus, eres = untw.pair_taus(ch, psi, sec)
            for c in CS:
                a = (phi[0] - c) / (phi[0] - phi[3]); kap = 1 - c * eps
                ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
                K = untw.Kval(ch, taus, kap, y)
                rhs = K * (untw.bra([a] * 3) @ psi) / (1 - kap * z[0]) ** 2
                lhs = untw.bra([ap, a, a]) @ psi
                out.append(dict(phi=phi, eps=eps, i=i, c=c, a=a, ap=ap, rhs=rhs, lhs=lhs, tgt=untw.C0[i][4], ov=ov))
        print('phi', phi, 'eps %.0e (%.0f s)' % (eps, time.time() - t0)); sys.stdout.flush()
pickle.dump(out, open('task3a.pkl', 'wb'))
for phi in PHIS:
    for c in CS:
        a = (phi[0] - c) / (phi[0] - phi[3])
        print('phi = %s, c = %.1f (a = %.4f):' % (phi, c, a))
        for i in STATES:
            rr = [r for r in out if r['phi'] == phi and r['c'] == c and r['i'] == i]
            print('   #%2d tgt=%+.6f%+.6fi  RHS(1e-4)=%+.6f%+.6fi  |RHS-tgt| = %s   |RHS-LHS| max %.0e   a\'-a(1e-4) = %.5f'
                  % (i, rr[0]['tgt'].real, rr[0]['tgt'].imag, rr[-1]['rhs'].real, rr[-1]['rhs'].imag,
                     ' '.join('%.1e' % abs(r['rhs'] - r['tgt']) for r in rr), max(abs(r['rhs'] - r['lhs']) for r in rr), (rr[-1]['ap'] - rr[-1]['a']).real))
