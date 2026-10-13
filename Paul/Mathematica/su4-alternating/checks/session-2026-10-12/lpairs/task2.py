# Task 2: even bridges at L = 3 through the untwisting limit, pair-1 word families + global torus (M = 3).
#  F1: kappa = 1/(z1 - a(z1-z4)) (a' - a = -1 exactly), a = 1/3   -> <W_a'(x)W_a^2|Psi> -> <W_-1 vac^2|Psi0> = C4 (M=1)
#  F2: kappa = 1/(z1 + a(z1-z4)) (a' - a = +1 exactly), a = 1/3   -> shift by -a': <vac W_-1^2|Psi0> = C2 (M=2)
#  F3: kappa = 1/z4, a = -1  (a' = 0)                               -> <vac W_-1^2|Psi0> = C2 (M=2)
#  F4: kappa = 1/z4, a = +1  (a' = 0) -> shift by -1                -> <W_-1 vac^2|Psi0> = C4 (M=1)
#  G : global rho(1-kappa Lambda), kappa as F1                      -> <W_-1^3|Psi0> = C0 (M=3)
#  formula: <W_a'(x)W_a^{L-1}|Psi> = K_kappa(tau) <W_a^L|Psi> / (1-kappa z1)^2
import numpy as np, sys, time, pickle
import untw
phi = [0.9, 0.35, 0.35, -1.1]; thdir = [0.0, 0.37, -0.52]
if len(sys.argv) > 1: phi = [float(x) for x in sys.argv[1].split(',')]
eps_list = (1e-2, 1e-3, 1e-4)
fams = {'F1': ('exact', -1, 1 / 3, 4), 'F2': ('exact', +1, 1 / 3, 2), 'F3': ('vac', None, -1.0, 2), 'F4': ('vac', None, 1.0, 4)}
rows = []
for eps in eps_list:
    t0 = time.time()
    ch, z = untw.make_chain(eps, phi, thdir)
    y = z[1]
    st = untw.twisted_states(ch)
    for i in sorted(st):
        psi, ov, sec, gap = st[i]
        taus, eres = untw.pair_taus(ch, psi, sec)
        rec = dict(eps=eps, i=i, dyn=tuple(untw.DYN[i]), E=untw.E0[i], ov=ov, gap=gap, eres=eres)
        for name, (mode, sgn, a, cidx) in fams.items():
            if mode == 'exact': kap = 1 / (z[0] - sgn * a * (z[0] - z[3]))
            else: kap = 1 / z[3]
            ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
            K = untw.Kval(ch, taus, kap, y)
            ovW = untw.bra([a] * 3) @ psi
            rhs = K * ovW / (1 - kap * z[0]) ** 2
            lhs = untw.bra([ap, a, a]) @ psi
            rec[name] = dict(kap=kap, ap=ap, K=K, ovW=ovW, rhs=rhs, lhs=lhs, tgt=untw.C0[i][cidx])
        # global, M = 3
        a = 1 / 3; kap = 1 / (z[0] + a * (z[0] - z[3])); ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
        Pg = ((1 - kap * z[3]) / (1 - kap * z[0])) ** 3 * untw.weight_factor(z, kap, untw.W0[i])
        ovW = untw.bra([a] * 3) @ psi
        rec['G'] = dict(kap=kap, ap=ap, K=Pg, ovW=ovW, rhs=Pg * ovW, lhs=untw.bra([ap] * 3) @ psi, tgt=untw.C0[i][0])
        rows.append(rec)
    print('eps = %.0e done (%.0f s)' % (eps, time.time() - t0)); sys.stdout.flush()
pickle.dump(rows, open('task2_%s.pkl' % ('_'.join('%g' % x for x in phi)), 'wb'))
for name in ('F1', 'F2', 'F3', 'F4', 'G'):
    print('==== family', name)
    for i in sorted(set(r['i'] for r in rows)):
        rr = [r for r in rows if r['i'] == i]
        tg = rr[0][name]['tgt']
        s = '  #%2d %s E=%7.3f tgt=%+.6f%+.6fi |' % (i, rr[0]['dyn'], rr[0]['E'], tg.real, tg.imag)
        for r in rr:
            d = r[name]
            s += ' eps=%.0e: ov=%.6f |RHS-LHS|=%.1e |RHS-tgt|=%.1e |<W|Psi>|=%.1e |K/N|=%.1e;' % (
                r['eps'], r['ov'], abs(d['rhs'] - d['lhs']), abs(d['rhs'] - d['tgt']), abs(d['ovW']), abs(d['K'] / (1 - d['kap'] * 1) ** 0 if name == 'G' else abs(d['K'])))
        print(s)
