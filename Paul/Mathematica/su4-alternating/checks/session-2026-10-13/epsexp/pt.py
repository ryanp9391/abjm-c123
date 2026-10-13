"""pt.py L [states] -- untwisted perturbation theory in eps for the twisted eigenvector Psi_eps = Psi0 + eps Psi1 + eps^2 Psi2 (weight sector of Psi0),
M(eps) = GEN-combination of t's (commuting family; same eigenvectors as every t_k(u)).  Computes
  omega_n = <vac|J^n|Psi_eps>/eps^n limit  (n = 1: <vac|J|Psi1> ; n = 2: <vac|J^2|Psi2>)
  X_{n-1}(u) = <vac|T14(u) J^{n-1}|Psi_{n-1}> and checks omega_n * (tau_vac - tau_Psi0)(u) = n (phi1-phi4) X_{n-1}(u)
  cluster (irrep) decomposition of the PT sum via the su(4) Casimir.
Also the first-order shifts of the pair eigenvalues X' by Hellmann-Feynman."""
import numpy as np, sys, math, pickle
from ue2 import *
np.set_printoptions(linewidth=220, precision=6, suppress=True)
Lp = int(sys.argv[1]); states = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
phi = [0.9, 0.35, 0.35, -1.1]; thdir = [0.0, 0.37, -0.52]
S = Setup(Lp); N = S.N; D = S.D
vac = S.vac()
def casimir(sec):
    C = None
    for i in range(4):
        for j in range(4):
            Eij = S.global_gen(i, j); Eji = S.global_gen(j, i)
            t = Eij @ Eji
            C = t if C is None else C + t
    return C
CAS = casimir(None)
sf = lambda x: complex(x)
ch0, z0 = S.chain(0.0, phi, thdir)
def Mof(eps, sec, ops=GEN):
    ch, z = S.chain(eps, phi, thdir)
    return S.sector_matrix(ch, ops, sec)
res = {}
which = states if states is not None else range(len(S.PSI0))
bysec = {}
for i in which: bysec.setdefault(tuple(S.W0[i]), []).append(i)
for w, ids in bysec.items():
    sec = S.SECT[w]; d = len(sec)
    M0 = Mof(0.0, sec); h1 = 1e-4; h2 = 1e-3
    V = (Mof(h1, sec) - Mof(-h1, sec)) / (2 * h1)
    M2 = (Mof(h2, sec) + Mof(-h2, sec) - 2 * M0) / (2 * h2 ** 2)
    ev, R = np.linalg.eig(M0); Rinv = np.linalg.inv(R)
    Csec = CAS[sec][:, sec].toarray() if hasattr(CAS[sec][:, sec], 'toarray') else CAS[sec][:, sec]
    # Casimir eigenvalue per eigenvector of M0 (valid when clusters are irreps): cas_m = (Rinv Csec R)_mm
    casm = np.real(np.diag(Rinv @ Csec @ R))
    # bras restricted to the sector
    vj = [(S.vacJn(k))[sec] for k in range(0, 4)]
    for i in ids:
        n = nval(Lp, S.DYN[i])
        psi0 = S.PSI0[i][sec]
        k0 = int(np.argmax(np.abs(Rinv @ psi0)))
        lam0 = ev[k0]; l0 = Rinv[k0]; l0 = l0 / (l0 @ psi0)
        # cluster size of lam0
        csize = int(np.sum(np.abs(ev - lam0) < 1e-8))
        Pm = np.eye(d, dtype=complex) - np.outer(psi0, l0) + np.outer(psi0, l0)  # identity
        Gres = np.linalg.inv(M0 - lam0 * np.eye(d) + np.outer(psi0, l0))
        lam1 = l0 @ (V @ psi0)
        psi1 = -Gres @ ((V - lam1 * np.eye(d)) @ psi0)
        rhs2 = (lam1 * np.eye(d) - V) @ psi1 - M2 @ psi0
        lam2 = -(l0 @ rhs2) * 0  # lam2 from solvability: l0.(rhs2 + lam2 psi0) = 0
        lam2 = -(l0 @ rhs2)
        psi2 = Gres @ (rhs2 + lam2 * psi0)
        out = dict(n=n, csize=csize, lam0=lam0, lam1=lam1, lam2=lam2)
        # vacuum J-overlaps order by order
        out['vacJ_psi'] = [[vj[k] @ p for p in (psi0, psi1, psi2)] for k in range(4)]
        # T14 bras: <vac|T14(u) J^k restricted to sector
        T14b = {}
        for u in (0.37, -0.29, 0.83):
            for k in range(0, 3):
                b = vac.reshape((4,) * N)
                bb = S.vacJn(k).reshape((4,) * N)
                T14b[(u, k)] = apply_t1G(ch0, u, bb, e41()).reshape(-1)[sec]
        out['T14b'] = T14b
        # untwisted tau of Psi0 and tau_vac at u
        tv = {}
        for u in (0.37, -0.29, 0.83):
            tv[u] = ((ch0.apply_t(1, u, vac.reshape((4,) * N)).reshape(-1) @ vac) / (vac @ vac), tau_ket(ch0, S.PSI0[i], 1, u))
        out['tv'] = tv
        # decomposition of X_{n-1}(u) over clusters (Casimir values)
        if n is not None and n >= 1:
            u = 0.37
            bra = T14b[(u, n - 1)]
            pk = [psi0, psi1][n - 1] if n <= 2 else None
            out['X'] = {u_: (T14b[(u_, n - 1)] @ (psi0 if n == 1 else psi1)) for u_ in (0.37, -0.29, 0.83)}
        # cluster decomposition of psi1 and of the contributions to omega_n
        comp1 = Rinv @ psi1
        contrib = {}
        for m in range(d):
            if abs(comp1[m]) > 1e-10 or True:
                key = (round(casm[m], 4))
                contrib.setdefault(key, 0)
        out['casm_psi1'] = {}
        for m in range(d):
            cm = round(float(casm[m]), 3)
            out['casm_psi1'].setdefault(cm, np.zeros(2, dtype=complex))
            if n is not None and n >= 1:
                bra = vj[n]            # <vac|J^n
                out['casm_psi1'][cm] += np.array([ (bra @ R[:, m]) * comp1[m] if n == 1 else (T14b[(0.37, n - 1)] @ R[:, m]) * comp1[m], 0])
        res[i] = out
        print('#%d %s n=%s dim=%d cluster(lam0)=%d  lam1=%s' % (i, S.DYN[i], n, d, csize, np.round(lam1, 6)))
        if n is not None and n >= 1:
            print('   <vac|J^k|Psi_m> (k=0..3 rows; m=0,1,2 cols):'); print(np.array(out['vacJ_psi']))
            for u in (0.37, -0.29, 0.83):
                tvu, tpu = tv[u]
                lhs = None
                print('   u=%.2f tau_vac-tau_Psi0 = %s ; X_{n-1}(u) = %s ; X/(tau_vac-tau)= %s' % (u, np.round(tvu - tpu, 6), np.round(out['X'][u], 6), np.round(out['X'][u] / (tvu - tpu), 6)))
            print('   Casimir clusters contributing to X_{n-1} = <vac|T14 J^{n-1}|Psi_{n-1}> via first-order PT (u=0.37):')
            for cm, v in sorted(out['casm_psi1'].items()):
                if abs(v[0]) > 1e-10: print('      Casimir %.3f : %s' % (cm, np.round(v[0], 6)))
pickle.dump(res, open('pt_L%d.pkl' % Lp, 'wb'))
