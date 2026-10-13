# generic theta, z2 = z3, L = 3, W_a refs: find a SPARSE localising combination producing the 5th localised direction
import numpy as np, itertools, time
from scipy.optimize import linprog
from pgbase import *
h = Fr(1,2)
zz = [Fr(3,2), Fr(5,4), Fr(5,4), Fr(4,7)]
thF = [Fr(-3,5), Fr(-11,8), Fr(-7,9)]; dA = [Fr(31,35), Fr(39,40), Fr(13,9)]
th = []
for f, d in zip(thF, dA): th += [f + d, f]
ch = Chain(['A','F']*3, th, h, diagm(zz))
def Wp(a): return [np.array([1,0,0,a], float), np.array([-a,0,0,1], float)]
avals = (1/3, -0.4, 1.7)
singles = [('F', k, m) for k in (1,2,3) for m in range(-2, 5)] + [('A', k, m) for k in (1,2,3) for m in range(-3, 2)]
def u_of(s): return (th[1] if s[0] == 'F' else th[0]) + s[2]*h
outs = {}
for r, a in enumerate(avals):
    ref = prod_state(Wp(a)*3, dtype=float); outs[(r, ())] = ref
    for s in singles: outs[(r, (s,))] = ch.apply_t(s[1], u_of(s), ref, reg=True).astype(float)
    for s1, s2 in itertools.combinations_with_replacement(singles, 2):
        outs[(r, (s1, s2))] = ch.apply_t(s2[1], u_of(s2), outs[(r, (s1,))], reg=True).astype(float)
words = [()] + [(s,) for s in singles] + list(itertools.combinations_with_replacement(singles, 2))
F = []; NL = []
for w in words:
    fp = []; npp = []
    for r, a in enumerate(avals):
        rest = prod_state(Wp(a)*2, dtype=float).reshape(-1)
        M = outs[(r, w)].reshape(16, -1); ps = M@rest/(rest@rest)
        fp.append(ps); npp.append((M - np.outer(ps, rest)).reshape(-1))
    F.append(np.concatenate(fp)); NL.append(np.concatenate(npp))
F = np.array(F).T; NL = np.array(NL).T
sc = np.maximum(np.abs(F).max(0), np.abs(NL).max(0)); sc[sc == 0] = 1
F = F/sc; NL = NL/sc
nz = np.linalg.norm(NL, axis=0) > 1e-9
indiv = [i for i in range(len(words)) if not nz[i]]
print('individually localising words:', [words[i] for i in indiv][:20], '... count', len(indiv))
B0 = F[:, indiv]; q, rr = np.linalg.qr(B0); rk0 = np.linalg.matrix_rank(B0, 1e-8)
print('rank of individually localised span:', rk0)
U, S, Vt = np.linalg.svd(NL, full_matrices=False)
rk = np.sum(S > 1e-9*S[0]); Aeq_nl = (S[:rk, None]*Vt[:rk])
null = Vt[rk:].T if False else None
# 5th direction: component of the localised image orthogonal to individual span
_, _, vh = np.linalg.svd(NL, full_matrices=True); nullsp = vh[rk:].T
img = F @ nullsp
Uq, Sq, _ = np.linalg.svd(B0, full_matrices=False); Q0 = Uq[:, Sq > 1e-8*Sq[0]]
perp = img - Q0 @ (Q0.T @ img)
up, sp_, _ = np.linalg.svd(perp, full_matrices=False)
print('dims of localised image outside the individual span:', np.sum(sp_ > 1e-8*sp_[0]), ' sv', sp_[:3])
d = up[:, 0]
# min L1:  NL c = 0,  P_perp F c = d   (c = cp - cm)
Pp = np.eye(F.shape[0]) - Q0 @ Q0.T
Aeq = np.vstack([Aeq_nl, Pp @ F]); beq = np.concatenate([np.zeros(rk), d])
n = F.shape[1]
res = linprog(np.ones(2*n), A_eq=np.hstack([Aeq, -Aeq]), b_eq=beq, bounds=(0, None), method='highs')
print(res.status, res.message)
c = res.x[:n] - res.x[n:]
supp = np.where(np.abs(c) > 1e-7*np.abs(c).max())[0]
print('support size', len(supp))
for i in supp: print('  ', words[i], '%.4e' % (c[i]/sc[i]))
print('check NL residual', np.linalg.norm(NL@c), ' perp component', np.linalg.norm(Pp@F@c - d))
