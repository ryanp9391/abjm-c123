# Task 3c: untwisting limit with z2 != z3: z = 1 + eps phi, then z2 += d/2, z3 -= d/2, d = eta*eps^q (q = 1, 2, 3).
# Word K_A (y -> (z2+z3)/2 linear, y^2 -> z2 z3).  Families F1 (a'-a = -1 exact, a = 1/3, tgt C4), F3 (kappa = 1/z4, a = -1, tgt C2).
import numpy as np, sys, time, pickle
import untw, eng, lp
h = 0.5; thdir = [0.0, 0.37, -0.52]; phi = (0.9, 0.35, 0.35, -1.1); eta = 0.5
STATES = [0, 4, 16, 17, 18, 20, 22, 25]
def KA(ch, taus, kap, z):
    p = 1; uF = ch.th[1]; u2 = uF + h
    l1, l2, l3 = eng.lam(ch, 1, uF, p), eng.lam(ch, 2, u2, p), eng.lam(ch, 3, u2, p)
    s = 1 - kap * (z[0] + z[3]); y1 = (z[1] + z[2]) / 2; y2 = z[1] * z[2]
    t1, t2, t3 = taus
    return s * s + kap * s * (t1 / (h**2 * l1) + t2 / (2 * h**4 * y1 * l2) + t3 / (6 * h**5 * y2 * l3)) + kap**2 * t3 * t1 / (24 * h**7 * y2 * l1 * l3)
out = []
for q in (1, 2, 3):
    for eps in (1e-2, 1e-3, 1e-4):
        t0 = time.time()
        def zfun(e, q=q):
            z = list(1 + e * np.array(phi)); d = eta * e ** q
            z[1] += d / 2; z[2] -= d / 2; return z
        ch, z = untw.make_chain(eps, phi, thdir, zfun=zfun)
        st = untw.twisted_states(ch, STATES)
        for i in STATES:
            psi, ov, sec, gap = st[i]
            taus, eres = untw.pair_taus(ch, psi, sec)
            for fam, a, cidx in (('F1', 1 / 3, 4), ('F3', -1.0, 2)):
                kap = 1 / (z[0] - a * (z[0] - z[3])) if fam == 'F1' else 1 / z[3]
                ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
                rhs = KA(ch, taus, kap, z) * (untw.bra([a] * 3) @ psi) / (1 - kap * z[0]) ** 2
                lhs = untw.bra([ap, a, a]) @ psi
                out.append(dict(q=q, eps=eps, i=i, fam=fam, rhs=rhs, lhs=lhs, tgt=untw.C0[i][cidx], ov=ov, gap=gap))
        print('q=%d eps=%.0e (%.0f s)' % (q, eps, time.time() - t0)); sys.stdout.flush()
pickle.dump(out, open('task3c.pkl', 'wb'))
for fam in ('F1', 'F3'):
    for q in (1, 2, 3):
        print('== %s, z2 - z3 = %.1f eps^%d' % (fam, eta, q))
        for i in STATES:
            rr = [r for r in out if r['fam'] == fam and r['q'] == q and r['i'] == i]
            print('   #%2d tgt=%+.6f%+.6fi | ov(1e-4)=%.6f | |LHS-tgt|: %s | |RHS-LHS|: %s | |RHS-tgt|: %s' % (
                i, rr[0]['tgt'].real, rr[0]['tgt'].imag, rr[-1]['ov'],
                ' '.join('%.1e' % abs(r['lhs'] - r['tgt']) for r in rr), ' '.join('%.1e' % abs(r['rhs'] - r['lhs']) for r in rr),
                ' '.join('%.1e' % abs(r['rhs'] - r['tgt']) for r in rr)))
