"""Physical/untwisting check: B(b) = <0| rho_12(1 + b E14^{diag}) ; C0 = lim <B(b)|Psi_eps> / (lam1^2 (-b)^(2-w1))."""
import numpy as np, sys
from core import *
from states import *
import wxf
D = wxf.load(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf')
def sparse_to_dense(sa):
    dims = int(sa.args[1][0]); rowptr, cols = sa.args[3][1]; vals = sa.args[3][2]
    v = np.zeros(dims, dtype=complex)
    for c, x in zip(np.array(cols).ravel(), vals): v[int(c) - 1] = x
    return v
L2 = D['L2']
data = [(np.array(s['w']), np.array(s['dyn']), s['E'], s['C'], sparse_to_dense(s['psi'])) for s in L2]
e1, e4 = np.eye(4)[0], np.eye(4)[3]
W = lambda a: np.kron(e1 + a * e4, e4 - a * e1)
vac = np.kron(e1, e4)
print('data check: <vac x W_-1|psi> vs C2, <W_-1 x vac|psi> (cyclic partner):')
for w, dyn, E, C, psi in data:
    print('  ', dyn, 'E=%.3f' % E, ' C2 =', np.round(C[2], 12), ' vac(x)W-1:', np.round(np.kron(vac, W(-1)) @ psi, 12), ' W-1(x)vac:', np.round(np.kron(W(-1), vac) @ psi, 12))
# random Borel n check (HW invariance) at the untwisted point
rng = np.random.default_rng(0)
n = np.eye(4) + np.triu(rng.normal(size=(4, 4)), 1)
def rho_bra(vs, g):  # A: v g ; F: u g^{-T}
    gi = np.linalg.inv(g).T
    return [v @ (g if types[s] == 'A' else gi) for s, v in enumerate(vs)]
B = rho_bra([e1 - e4, e4 + e1, e1, e4], n)
print('HW invariance, random n: ', [np.round(kronvecs(Ring('float'), B) @ d[4], 12) for d in data])
phi = [Fr(9, 10), Fr(2, 5), Fr(-3, 10), Fr(-11, 10)]
thA = [Fr(1), Fr(0), Fr(1), Fr(0)]
results = []
for eps in (Fr(1, 100), Fr(1, 1000), Fr(1, 10000)):
    dl = eps
    th = [thA[s] + dl * th0[s] for s in range(4)]
    z = [1 + eps * phi[j] for j in range(4)]
    Rf = Ring('float')
    chd = Chain(Rf, th=th, G=[[z[i] if i == j else 0 for j in range(4)] for i in range(4)])
    evs, Vd = eigenstates_diag(chd)
    Wt = weights()
    lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
    lamf = np.array([float(x) for x in lam]); one = np.ones(4)
    for b in (0.7, -1.3):
        Bb = np.kron(np.kron(lamf + b * lamf[0] * e4, one - b * e1), np.kron(lamf, one))
        for k, (w, dyn, E, C, psi) in enumerate(data):
            sec = [i for i in range(256) if tuple(Wt[i]) == tuple(w)]
            # candidate eigenstates: columns of Vd supported in sector
            cand = [c for c in range(256) if np.linalg.norm(Vd[sec, c]) > 0.999 * np.linalg.norm(Vd[:, c])]
            ov = [abs(np.vdot(psi, Vd[:, c])) / np.linalg.norm(Vd[:, c]) for c in cand]
            c = cand[int(np.argmax(ov))]
            Psi = Vd[:, c] / np.linalg.norm(Vd[:, c])
            ph = np.vdot(psi, Psi); Psi = Psi * abs(ph) / ph
            w1 = int(w[0])
            Cpred = (Bb @ Psi) / (lamf[0] ** 2 * (-b) ** (2 - w1))
            results.append((float(eps), b, k, max(ov), Cpred, C[2]))
            print(f'eps=delta={float(eps):.0e} b={b:+.1f} state {k} dyn={dyn} |<psi0|Psi>|={max(ov):.6f}  C_pred={np.round(Cpred, 8)}  data C2={np.round(C[2], 8)}')
