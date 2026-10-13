"""phi- and y0-dependence of the untwisting limit of the A-site (odd-bridge) word on <0|: [1,0,1] state, L = 2."""
import numpy as np
from untw import *
L = 2; S = Setup(L); DATA = load_data_states(L)
st = DATA[2]; vacst = DATA[0]
def run(phi, thdir, eps, y0):
    ch, th, z = S.chain(eps, phi, thdir)
    lam = lam_of(z); c4 = np.prod(z); u0 = th[1] + h
    T3 = S.op(ch, 3, u0, reg=True)
    Mg = S.op(ch, 1, th[1]) + 0.61 * S.op(ch, 3, -0.29) + 0.23 * S.op(ch, 2, 0.83) + 0.3 * S.op(ch, 1, 0.71)
    bra0 = prod_state([lam, np.ones(4)] * L, dtype=complex).reshape(-1)
    B0T3 = bra0 @ T3
    # companion e4 coefficient of <0|X3 on pair 1: read off from vacuum overlap structure: <0|X3 = c3 * (e4 e4)(x) vac
    Vi_T, V = comp_to_diag_site(z)
    pred = np.kron(np.kron(np.kron(np.eye(4)[3] @ Vi_T, np.eye(4)[3] @ V), lam), np.ones(4))
    c3 = (B0T3 @ pred.conj()) / (pred @ pred.conj())
    alpha = c4 / z[0] + 1 / (lam[0] * y0)
    B = alpha * bra0 + B0T3 / c3
    psi, psiL, ov, gap = S.twisted_pair(Mg, st['w'], st['psi'])
    return (B @ psi) / (B @ vacst['psi'])
for phi in ([0.9, 0.35, -0.3, -1.1], [0.2, -0.7, 0.5, 1.3], [1.0, 0.0, -0.5, 0.25]):
    for y0 in (0.7, -1.3):
        vals = [run(phi, [0.0, 0.37], eps, y0) for eps in (1e-3, 1e-4)]
        print('phi', phi, 'y0', y0, ' C(eps=1e-3,1e-4) =', np.round(vals, 6), ' |C|/|c_y y0| =', np.round(np.abs(vals) / (0.288675134594813 * abs(y0)), 6))
