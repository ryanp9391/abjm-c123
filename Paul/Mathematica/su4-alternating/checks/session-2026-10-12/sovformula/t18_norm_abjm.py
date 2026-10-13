"""ABJM point (pair by pair), L = 2, companion twist: complete regularised bra and ket SoV bases, Gram matrix, measure,
and the SoV norm 1/N = sum_{x,y} tau~_y M_yx tau_x checked on joint eigenstates (eigenvalue products vs direct)."""
import time, numpy as np
from regbases import *
from collections import Counter
T0 = time.time()
R = Reg()
Wb = R.bra_words(); Yk = R.ket_words()
print('ranks bra/ket (mod p):', S.rank_mod_p(Wb), S.rank_mod_p(Yk.T), ' %.0f s' % (time.time() - T0), flush=True)
G = S.gram(Wb, Yk)
nz = [(S.WORDS[i], S.WORDS[j]) for i in range(256) for j in range(256) if G[i, j] != 0]
print('Gram nonzeros:', len(nz), ' charge sums:', sorted(Counter(sum(x) + sum(y) for x, y in nz).items()), flush=True)
print('Gram triangular (x <-> 3-x relabelled, acyclic):', S.is_triangular(G))
Gf = np.array(G, dtype=float); Mf = np.linalg.inv(Gf)        # M[y, x] with G[x, y] = <x|y>
print('measure nonzeros (|M| > 1e-12 max):', int(np.sum(np.abs(Mf) > 1e-12 * np.abs(Mf).max())), ' cond(G) = %.2e' % np.linalg.cond(Gf))
# joint eigenstates (companion frame): diagonalise a generic combination (bra action matrices; kets: M psi)
th = R.th; h = R.h
Hc = sum(c * np.array(R.op(a, u), dtype=float) for a, u, c in ((1, Fr(37, 100), 1.0), (2, Fr(-113, 100), -0.7), (3, Fr(3, 2), 0.3), (1, Fr(-71, 100), 1.1)))
ev, VR = np.linalg.eig(Hc); VL = np.linalg.inv(VR)
vac = np.array(S.sov_vacuum(), dtype=float); ket0 = np.array(S.product_vec([S.unit(i) for i in S.KET_VAC]), dtype=float)
Wf = np.array(Wb, dtype=float); Yf = np.array(Yk, dtype=float)
errs = []; gaps = np.sort(np.abs(ev[:, None] - ev[None, :]) + np.eye(256) * 1e9, axis=1)[:, 0]
for n in range(0, 256, 9):
    psi = VR[:, n]; psiL = VL[n]
    tau = (Wf @ psi) / (vac @ psi)              # <x|Psi>/<0|Psi>
    taut = (psiL @ Yf) / (psiL @ ket0)           # <PsiL|y>/<PsiL|0bar>
    invN = taut @ Mf @ tau
    direct = (psiL @ psi) / ((vac @ psi) * (psiL @ ket0))
    errs.append(abs(invN / direct - 1))
print('1/N = sum tau~ M tau vs <PsiL|Psi>/(<0|Psi><PsiL|0bar>): max rel err %.1e over %d states (min gap %.1e)' % (max(errs), len(errs), gaps.min()))
# eigenvalue-product form of tau_x (bra) for a few states: pair factors lambda(a,k) tauF^k
for n in (0, 77, 151):
    psi = VR[:, n]
    def eig(M): v = np.array(M, dtype=float) @ psi; return (psi.conj() @ v) / (psi.conj() @ psi)
    P = [R.pair_ops_bra(0), R.pair_ops_bra(1)]
    lam = []
    for p in (0, 1):
        tF = eig(P[p]['t1F']); x = eig(P[p]['O'][1]); w = eig(P[p]['O'][2]); y = eig(P[p]['O'][3]); t2 = eig(P[p]['t2u0'])
        lam.append((tF, [1, x, w, y], t2))
    tau = (Wf @ psi) / (vac @ psi)
    e = 0
    for i, (a1, k2, a3, k4) in enumerate(S.WORDS):
        f = 1
        for p, (a, k) in enumerate(((a1, k2), (a3, k4))):
            tF, L4, t2 = lam[p]
            f *= (t2 if (a == 2 and k == 3) else L4[a]) * tF ** k
        e = max(e, abs(f - tau[i]) / max(abs(tau[i]), 1e-300))
    print('state', n, ': tau_x = prod_pairs lambda(a,k) tauF^k, max rel err %.1e' % e)
print('time %.0f s' % (time.time() - T0))
import pickle; pickle.dump(dict(G=G), open('gram_abjm_L2.pkl', 'wb'))
