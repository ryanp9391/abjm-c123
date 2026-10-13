import numpy as np
from untw import *
L = 2; S = Setup(L); DATA = load_data_states(L)
st = DATA[2]; vacst = DATA[0]
phi = [0.9, 0.35, -0.3, -1.1]; thdir = [0.0, 0.37]; y0 = 0.7
for eps in (1e-2, 1e-3):
    ch, th, z = S.chain(eps, phi, thdir)
    lam = lam_of(z); c4 = np.prod(z); u0 = th[1] + h
    T3 = S.op(ch, 3, u0, reg=True)
    Mg = S.op(ch, 1, th[1]) + 0.61 * S.op(ch, 3, -0.29) + 0.23 * S.op(ch, 2, 0.83) + 0.3 * S.op(ch, 1, 0.71)
    bra0 = prod_state([lam, np.ones(4)] * L, dtype=complex).reshape(-1)
    B0T3 = bra0 @ T3
    Vi_T, V = comp_to_diag_site(z)
    pred = np.kron(np.kron(np.kron(np.eye(4)[3] @ Vi_T, np.eye(4)[3] @ V), lam), np.ones(4))
    c3 = (B0T3 @ pred.conj()) / (pred @ pred.conj())
    l3 = 1.0
    for s in range(2, 4):
        t = th[s]; k = 3
        fac = [u0 - t - j * h for j in range(k)] if s % 2 == 0 else [u0 - t - k * h] + [u0 - t - j * h for j in range(k - 1)]
        for f in fac:
            if abs(f) > 1e-14: l3 *= f
    print('eps', eps, 'c3 =', c3, ' -12h^5 l3 =', -12 * h ** 5 * l3, ' resid', np.linalg.norm(B0T3 - c3 * pred) / np.linalg.norm(B0T3))
    alpha = c4 / z[0] + 1 / (lam[0] * y0)
    psi, psiL, ov, gap = S.twisted_pair(Mg, st['w'], st['psi'])
    psv, _, _, _ = S.twisted_pair(Mg, vacst['w'], vacst['psi'])
    for name, B in (('fit', alpha * bra0 + B0T3 / c3), ('formula', alpha * bra0 - B0T3 / (12 * h ** 5 * l3))):
        print('   ', name, ' <B|Psi> =', B @ psi, ' <B|vac_tw> =', B @ psv, ' <B|vac_data> =', B @ vacst['psi'], ' ratio', (B @ psi) / (B @ psv))
    # pure components
    print('    <0|Psi> =', bra0 @ psi, '  <0|X3|Psi> =', B0T3 @ psi, ' alpha =', alpha, ' |Psi_vac_tw - data|', np.linalg.norm(psv - vacst['psi']))
