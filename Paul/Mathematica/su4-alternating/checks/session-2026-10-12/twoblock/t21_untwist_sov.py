"""SoV side of the untwisting check: sum_x c_x(b) W_x <0|Psi> (companion frame) vs direct overlap and data."""
import numpy as np
from core import *
from states import *
exec(open('t20_untwist.py').read().split("phi = [")[0].split("print('data check")[0])
phi = [Fr(9, 10), Fr(2, 5), Fr(-3, 10), Fr(-11, 10)]
thA = [Fr(1), Fr(0), Fr(1), Fr(0)]
for eps in (Fr(1, 100), Fr(1, 1000)):
    th = [thA[s] + eps * th0[s] for s in range(4)]
    z = [1 + eps * phi[j] for j in range(4)]
    Rf = Ring('float')
    chd = Chain(Rf, th=th, G=[[z[i] if i == j else 0 for j in range(4)] for i in range(4)])
    chc = Chain(Rf, th=th, z=z)
    evs, Vd = eigenstates_diag(chd)
    Rm = Rmat(z); V = Rm @ Vd
    M, vac0 = sov_words(chc)
    print('eps=delta=%g  cond(M)=%.2e' % (float(eps), np.linalg.cond(M)))
    Ws = wavevals(chc, V); Wx = Wprod(Ws)
    lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
    lamf = np.array([float(x) for x in lam]); one = np.ones(4)
    Vf = np.array([[float(x) for x in row] for row in Vand(z)]); Vi = np.linalg.inv(Vf)
    vz4 = Vf[:, 3]; l1 = Vi[0]
    Wt = weights()
    for b in (0.7,):
        Bc = np.kron(np.kron(e1 + b * lamf[0] * vz4, e4 - b * l1), np.kron(e1, e4))   # companion frame
        Bd = np.kron(np.kron(lamf + b * lamf[0] * e4, one - b * e1), np.kron(lamf, one))  # diag frame
        print('   frame map check |Bc R - Bd|/|Bd| =', np.linalg.norm(Bc @ Rm - Bd) / np.linalg.norm(Bd))
        c = np.linalg.solve(M.T, Bc)
        print('   nnz(|c|>1e-9 max):', int(np.sum(np.abs(c) > 1e-9 * np.abs(c).max())), ' residual |cM-B|/|B|', np.linalg.norm(c @ M - Bc) / np.linalg.norm(Bc))
        for k, (w, dyn, E, C, psi) in enumerate(data):
            sec = [i for i in range(256) if tuple(Wt[i]) == tuple(w)]
            cand = [cc for cc in range(256) if np.linalg.norm(Vd[sec, cc]) > 0.999 * np.linalg.norm(Vd[:, cc])]
            ov = [abs(np.vdot(psi, Vd[:, cc])) / np.linalg.norm(Vd[:, cc]) for cc in cand]
            cc = cand[int(np.argmax(ov))]
            nrm = np.linalg.norm(Vd[:, cc]); ph = np.vdot(psi, Vd[:, cc]); ph = ph / abs(ph)
            direct = (Bd @ Vd[:, cc]) / nrm / ph
            sov = (c @ Wx[:, cc]) * (vac0 @ V[:, cc]) / nrm / ph
            w1 = int(w[0]); norm = lamf[0] ** 2 * (-b) ** (2 - w1)
            print(f'   state {k} {dyn}: direct/norm={np.round(direct / norm, 8)}  SoV/norm={np.round(sov / norm, 8)}  rel.diff={abs(sov - direct) / max(abs(direct), 1e-300):.1e}  data C2={np.round(C[2], 8)}  <0|Psi>/lam1^2={abs(vac0 @ V[:, cc] / nrm) / lamf[0] ** 2:.2e}')
