import pickle, numpy as np, sys, math
from ue_base import *
Lp = int(sys.argv[1]); d = pickle.load(open('scan2_L%d.pkl' % Lp, 'rb')); rows = d['rows']; a = d['a']; c = d['c']; phi = d['phi']
S = Setup(Lp)
np.set_printoptions(linewidth=220, precision=6, suppress=True)
def nval_(i): return nval(Lp, S.DYN[i])
print('L=%d a=%.4f c=%.2f phi=%s' % (Lp, a, c, phi))
for i in sorted(rows):
    rr = rows[i]; n = nval_(i)
    eps = np.array([r['eps'] for r in rr])
    print('---- state #%d %s n=%s E=%.3f' % (i, S.DYN[i], n, S.E0[i]))
    # identity check
    mx = 0
    for r in rr:
        for (u, tv, tp, T14, T14v) in r['ident']:
            lhs = (tv - tp) * r['w']; rhs = a * (0.9 - (-1.1)) * 0  # placeholder
    # fit sign/value: lhs / (a eps (phi1-phi4) T14)
    r = rr[-1]; eps_ = r['eps']
    for (u, tv, tp, T14, T14v) in r['ident']:
        lhs = (tv - tp) * r['w']
        z1mz4 = eps_ * (phi[0] - phi[3])
        print('   u=%.2f: (tvac - tau) w = %+.3e%+.3ei ; (z1-z4)<W|T14|Psi> = %+.3e%+.3ei ; ratio %+.6f%+.6fi (a=%.4f)' % (u, lhs.real, lhs.imag, (z1mz4 * T14).real, (z1mz4 * T14).imag, *( (lhs / (z1mz4 * T14)).real, (lhs/(z1mz4*T14)).imag ), a) if abs(T14) > 1e-14 else '   u=%.2f: T14 = 0; lhs = %.2e' % (u, abs(lhs)))
