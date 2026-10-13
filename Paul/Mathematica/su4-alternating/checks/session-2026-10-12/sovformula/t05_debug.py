import numpy as np
from untw import *
L = 2; S = Setup(L); DATA = load_data_states(L)
idx = lambda a, b, c, d: np.ravel_multi_index((a - 1, b - 1, c - 1, d - 1), (4,) * 4)
for st in DATA:
    p = st['psi']
    print(st['dyn'], st['w'], ' <e1e4e4e4|Psi0> =', np.round(p[idx(1, 4, 4, 4)], 6), ' <e4e4e1e4|Psi0> =', np.round(p[idx(4, 4, 1, 4)], 6),
          ' <e1e4e1e4|> =', np.round(p[idx(1, 4, 1, 4)], 6))
phi = [0.9, 0.35, -0.3, -1.1]; thdir = [0.0, 0.37]
eps = 1e-3
ch, th, z = S.chain(eps, phi, thdir)
lam = lam_of(z); c4 = np.prod(z); u0 = th[1] + h
T3 = S.op(ch, 3, u0, reg=True)
Mg = S.op(ch, 1, th[1]) + 0.61 * S.op(ch, 3, -0.29) + 0.23 * S.op(ch, 2, 0.83) + 0.3 * S.op(ch, 1, 0.71)
bra0 = prod_state([lam, np.ones(4)] * L, dtype=complex).reshape(-1)
B3 = bra0 @ T3
# expected: companion -12 h^5 lam3 e4 e4 (x) e1 e4 -> diag
Vi_T, V = comp_to_diag_site(z)
pred = np.kron(np.kron(np.kron(np.eye(4)[3] @ Vi_T, np.eye(4)[3] @ V), np.eye(4)[0] @ Vi_T), np.eye(4)[3] @ V)
c = (B3 @ pred.conj()) / (pred @ pred.conj())
print('B3 = c * pred, c =', c, ' residual', np.linalg.norm(B3 - c * pred) / np.linalg.norm(B3))
for st in DATA[2:3]:
    psi, psiL, ov, gap = S.twisted_pair(Mg, st['w'], st['psi'])
    print('ov', ov, 'gap', gap)
    print('  twisted <e4e4e1e4|Psi>', psi[idx(4, 4, 1, 4)], ' <e1e4e4e4|Psi>', psi[idx(1, 4, 4, 4)], ' psi0:', st['psi'][idx(4, 4, 1, 4)])
    n, Dp = nprime(z)
    g = n @ Dp
    # physical bra of pred: pair vectors times g^{-1} (A), g^T (F)
    gi = np.linalg.inv(g)
    vA = np.eye(4)[3] @ Vi_T @ gi; uF = np.eye(4)[3] @ V @ g.T
    print('  physical A (from companion e4):', np.round(vA, 6), ' F:', np.round(uF, 6))
    print('  <pred|Psi> =', pred @ psi, '  pred_phys . psi0 (no Dw)=', np.kron(np.kron(np.kron(vA, uF), np.eye(4)[0] @ Vi_T @ gi), np.eye(4)[3] @ V @ g.T) @ st['psi'])
