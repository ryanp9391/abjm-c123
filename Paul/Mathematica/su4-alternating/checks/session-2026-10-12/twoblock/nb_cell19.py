import numpy as np
from core import *
from nb_common import *
import wxf
D = wxf.load(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf')
def sp2d(sa):
    dims = int(sa.args[1][0]); cols = sa.args[3][1][1]; vals = sa.args[3][2]
    v = np.zeros(dims, dtype=complex)
    for c, x in zip(np.array(cols).ravel(), vals): v[int(c) - 1] = x
    return v
data = [(tuple(int(x) for x in s['w']), [int(x) for x in s['dyn']], s['E'], s['C'], sp2d(s['psi'])) for s in D['L2']]
del D
phi = [Fr(9, 10), Fr(2, 5), Fr(-3, 10), Fr(-11, 10)]; thA = [Fr(1), Fr(0), Fr(1), Fr(0)]
e1, e4 = np.eye(4)[0], np.eye(4)[3]
def match(Vd, w, psi):
    sec = [i for i in range(256) if WT[i] == w]
    cand = [c for c in range(256) if np.linalg.norm(Vd[sec, c]) > (1 - 1e-6) * np.linalg.norm(Vd[:, c])]
    ov = [abs(np.vdot(psi, Vd[:, c])) / np.linalg.norm(Vd[:, c]) for c in cand]
    c = cand[int(np.argmax(ov))]
    f = 1 / np.linalg.norm(Vd[:, c]); ph = np.vdot(psi, Vd[:, c]); f *= abs(ph) / ph
    return c, max(ov), f
print('| ε = δ | b | Dynkin | E_MZ | \|⟨ψ₀\|Ψ_ε⟩\| | C_pred | data C₂ | \|C_pred − C₂\| |'); print('|---|---|---|---|---|---|---|---|')
for eps in (Fr(1, 100), Fr(1, 1000), Fr(1, 10000)):
    th = [thA[s] + eps * th0[s] for s in range(4)]; z = [1 + eps * phi[j] for j in range(4)]
    chd = Chain(Ring('float'), th=th, G=[[z[i] if i == j else 0 for j in range(4)] for i in range(4)])
    Vd, gap = states_diag(chd)
    lam = np.array([float(Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j])) for j in range(4)]); one = np.ones(4)
    for (w, dyn, E, C, psi) in data:
        c, ov, f = match(Vd, w, psi)
        for b in (0.7, -1.3):
            Bd = np.kron(np.kron(lam + b * lam[0] * e4, one - b * e1), np.kron(lam, one))
            cp = f * (Bd @ Vd[:, c]) / (lam[0] ** 2 * (-b) ** (2 - w[0]))
            cd = C[2] if not isinstance(C[2], complex) else C[2]
            cps = ('%.8f' % cp.real) + ('' if abs(cp.imag) < 1e-12 else ' %+.8fi' % cp.imag)
            print(f'| {float(eps):.0e} | {b:+.1f} | {dyn} | {E:.4f} | {ov:.8f} | {cps} | {complex(cd).real:.8f} | {abs(cp - cd):.2e} |')
