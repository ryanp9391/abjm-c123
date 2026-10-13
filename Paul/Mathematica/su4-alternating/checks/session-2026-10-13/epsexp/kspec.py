"""kspec.py -- spectrum of the untwisted K_1 operator (z = 1, kappa = 1, y = 1, theta_F = 0, theta_A = 2h) on whole weight sectors (ALL eigenstates,
not only the HW ones): K_1 = 1 + c1 X1 + c2 X2 + c3 X3 + c4 X3 X1 (Kcoef at kappa = 1, s = -1) with X1 = t1(theta_F), X2 = t2(u0), X3 = t3hat(u0), u0 = theta_F + h.
Also the su(4) Casimir per eigenvector to name the irrep, and A+B+C, A*C (A = tauF/(h^2 l1) ...)."""
import numpy as np, sys
from ue2 import *
np.set_printoptions(linewidth=220, precision=5, suppress=True)
Lp = int(sys.argv[1]); secs = [tuple(int(x) for x in s.split(',')) for s in sys.argv[2:]]
S = Setup(Lp); ch0 = untwisted_chain(Lp)
cc = S.Kcoef(ch0, 1.0, 1.0)
uF = ch0.th[1]; u0 = uF + h
def Kmat(sec):
    M1 = S.sector_matrix(ch0, [(1.0, 1, uF, False)], sec); M2 = S.sector_matrix(ch0, [(1.0, 2, u0, False)], sec); M3 = S.sector_matrix(ch0, [(1.0, 3, u0, True)], sec)
    d = len(sec); I = np.eye(d)
    return cc[0] * I + cc[1] * M1 + cc[2] * M2 + cc[3] * M3 + cc[4] * (M3 @ M1), M1, M3
CAS = None
for w in secs:
    sec = S.SECT[w]; K, M1, M3 = Kmat(sec); d = len(sec)
    ev, R = np.linalg.eig(K)
    # semisimplicity check: condition of R
    cnd = np.linalg.cond(R)
    # commutation with Casimir to label irreps
    if CAS is None:
        CAS = None
        for i in range(4):
            for j in range(4):
                t = S.global_gen(i, j) @ S.global_gen(j, i)
                CAS = t if CAS is None else CAS + t
    Cs = CAS[sec][:, sec].toarray()
    Rinv = np.linalg.inv(R); cas = np.real(np.diag(Rinv @ Cs @ R))
    # also tau_F eigenvalue per eigenvector (momentum)
    uF_ev = np.diag(Rinv @ M1 @ R)
    nz = np.sum(np.abs(ev) < 1e-7)
    print('sector %s  L=%d dim=%d  cond(R)=%.1e  #(K1=0) = %d' % (w, Lp, d, cnd, nz))
    vals = {}
    for e, c_ in zip(ev, cas):
        key = (round(float(c_), 2), round(e.real, 4), round(e.imag, 4))
        vals[key] = vals.get(key, 0) + 1
    for k_, v in sorted(vals.items()): print('    Casimir %.2f  K1 = %+.4f%+.4fi   multiplicity %d' % (k_[0], k_[1], k_[2], v))
