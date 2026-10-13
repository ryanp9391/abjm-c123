"""Task 3 (L = 2, numerics): SoV-normalised formula through the untwisting limit.
Reference <0|_comp = (e1 e4)^L, diag frame (lambda (x) 1)^L; ket reference |0bar>_diag = (1 (x) lambda)^L.
For each HW state: twisted right/left eigenvectors, SoV norm  N = <0|Psi><PsiL|0bar>/<PsiL|Psi>,  weight factor lambda^w,
check |<0|Psi>|^2/||Psi||^2 = |lambda^w N|  (PsiL -> Psi^dagger), and the A-site (odd-bridge) word K = alpha + beta X3:
   <0|K|Psi> = P(tau) <0|Psi>,  P = alpha + beta*tau3hat(thF+h),
   C_formula = [P(tau)/P(tau_vac)] * sqrt(lambda^(w - w_vac) N / N_vac)   vs  physical F/F_vac (twoblock forms)."""
import numpy as np, sys, json
from untw import *
L = 2
S = Setup(L)
DATA = load_data_states(L)
print('HW states (python):', len(S.hw), ' data states:', len(DATA))
forms = json.load(open(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-twoblock-forms.json'))['L2']
phi = [0.9, 0.35, -0.3, -1.1]; thdir = [0.0, 0.37]
y0 = 0.7            # target physical y = p32/(a3 d2) of the rotated pair
def lam_hat(th, k, u0):
    r = 1.0
    for s in range(2, 2 * L):
        t = th[s]
        fac = [u0 - t - j * h for j in range(k)] if s % 2 == 0 else [u0 - t - k * h] + [u0 - t - j * h for j in range(k - 1)]
        for f in fac:
            if abs(f) > 1e-14: r *= f
    return r
for eps in (1e-2, 1e-3, 1e-4):
    ch, th, z = S.chain(eps, phi, thdir)
    lam = lam_of(z); c4 = np.prod(z)
    uF, u0 = th[1], th[1] + h
    T1 = S.op(ch, 1, uF); T2 = S.op(ch, 2, u0); T3 = S.op(ch, 3, u0, reg=True)
    Mg = T1 + 0.61 * S.op(ch, 3, -0.29) + 0.23 * S.op(ch, 2, 0.83) + 0.3 * S.op(ch, 1, 0.71)
    bra0 = prod_state([lam, np.ones(4)] * L, dtype=complex).reshape(-1)
    ket0b = prod_state([np.ones(4), lam] * L, dtype=complex).reshape(-1)
    # A-site word: companion pair-1 A vector e1 -> alpha e1 + b' e4 with alpha/b' = c4/z1 + 1/(lam1 y0)  (y = y0 exactly)
    l3 = lam_hat(th, 3, u0)
    bp = 1.0; alpha = bp * (c4 / z[0] + 1 / (lam[0] * y0))
    beta = -bp / (12 * h ** 5 * l3)          # <0|X3 = -12 h^5 lam3 e4 (x) e4 (x) vac  (companion)
    B = alpha * bra0 + beta * (bra0 @ T3)
    rows = []
    for st, fm in zip(DATA, forms):
        w, psi0 = st['w'], st['psi']
        psi, psiL, ov, gap = S.twisted_pair(Mg, w, psi0)
        t3 = (psi @ T3.T @ psi.conj()) / (psi @ psi.conj()) if False else None
        # eigenvalue of X3 from the ket action: (T3^T psi) = tau psi
        v = T3 @ psi; tau3 = np.vdot(psi, v) / np.vdot(psi, psi)
        P = alpha + beta * tau3
        o0 = bra0 @ psi; oB = B @ psi; o0b = psiL @ ket0b; nn = psiL @ psi
        Nsov = o0 * o0b / nn
        lw = np.prod([lam[j] ** w[j] for j in range(4)])
        rows.append(dict(w=w, ov=ov, gap=gap, P=P, o0=o0, oB=oB, Nsov=Nsov, lw=lw, res=abs(oB - P * o0) / max(abs(oB), 1e-300),
                         herm=abs(o0) ** 2 / abs(lw * Nsov), dyn=fm['dyn']))
    vac = [r for r in rows if r['w'] == (2, 0, 0, -2)][0]
    print(f'eps = {eps:.0e}')
    for r in rows:
        Dw = (lam[1] / lam[0]) ** r['w'][1] * (lam[2] / lam[0]) ** r['w'][2]
        Cf = (r['P'] / vac['P']) * np.sqrt(r['lw'] / vac['lw'] * r['Nsov'] / vac['Nsov'] + 0j) / Dw
        Cdir = r['oB'] / vac['oB'] / Dw
        print(f"  dyn={r['dyn']} w={r['w']} ov={r['ov']:.6f} |<0|Psi>|^2/|lam^w N|={r['herm']:.8f}  "
              f"<B|Psi>=P<0|Psi> res {r['res']:.1e}  |C_formula|={abs(Cf):.8f}  |<B|Psi>/<B|vac>|={abs(Cdir):.8f}")
    print('   physical prediction |c_y y0| for [1,0,1] = %.8f, others 0, vacuum 1' % (0.288675134594813 * y0))
