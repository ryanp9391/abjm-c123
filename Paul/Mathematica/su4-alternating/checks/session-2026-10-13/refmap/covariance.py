"""covariance.py -- item (6): GL(4) covariance <vac|rho(g) X[t^Lambda] = <vac| X[t^{G'}] rho(g); admissibility in terms of G'.
Also item (d): Lie algebra of pair-1 rotations in the 5-word span.  Python fallback for notebook sections 6, 9.  Output: covariance_out.txt"""
from rm_base import *
from rm_base import _OPS
from helpers import pred, tw_u
import itertools
out = open('covariance_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
y = Fr(5, 4); zz = [Fr(3, 2), y, y, Fr(4, 7)]
Lam = diagm(zz).tolist() if hasattr(diagm(zz), 'tolist') else diagm(zz)
Lam = [[Fr(x) for x in r] for r in Lam]

def rhoM(g):
    """rho(g) = g (x) g^{-T} as 16x16 (bra action psi -> psi.rho)"""
    gi = tr(inv(g))
    return [[g[i][k] * gi[j][l] for k in range(4) for l in range(4)] for i in range(4) for j in range(4)]

# (a) which G' makes rho(g) t_k^Lambda rho(g)^-1 = t_k^{G'} ?  (k = 1,2,3, one pair, generic z, theta = (1,0), u = 37/100)
zg = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
g = [[Fr(2), Fr(1, 3), Fr(-1), Fr(2)], [Fr(0), Fr(1), Fr(3, 2), Fr(1, 5)], [Fr(0), Fr(0), Fr(5, 3), Fr(-2)], [Fr(0), Fr(0), Fr(0), Fr(3, 4)]]
Lg = [[Fr(x) if i == j else Fr(0) for j, x in enumerate(zg)] for i in range(4)]
gi = inv(g)
cands = {'g^-1 L g': mm(mm(gi, Lg), g), 'g L g^-1': mm(mm(g, Lg), gi), 'g^T L g^-T': mm(mm(tr(g), Lg), tr(gi)), 'g^-T L g^T': mm(mm(tr(gi), Lg), tr(g))}
from pgbase import Chain
def ops_for(G):
    ch = Chain(['A', 'F'], [Fr(1), Fr(0)], h, np.array(G, dtype=object))
    res = []
    for k, u, reg in [(1, Fr(37, 100), False), (2, Fr(37, 100), False), (3, Fr(37, 100), False)]:
        M = []
        for a in range(4):
            for b in range(4):
                psi = np.zeros((4, 4), dtype=object); psi[:] = Fr(0); psi[a, b] = Fr(1)
                M.append([Fr(x) for x in ch.apply_t(k, u, psi, reg=reg).reshape(16)])
        res.append(M)
    return res
Rg = rhoM(g); Rgi = inv(Rg)
base = ops_for(Lg)
for name, G in cands.items():
    oc = ops_for(G)
    P('(6a) rho(g) t_k^Lambda rho(g)^-1 == t_k^{G\'} with G\' =', name, ':', [mm(mm(Rg, base[k]), Rgi) == oc[k] for k in range(3)])

# (b) the covariance of K and the statement <vac|K[G'] ~ <vac|rho(1-kappa G') for families
def Kop(kap, G, zz):
    ops = pair_ops(np.array(G, dtype=object))
    yy = zz[1]; s = 1 - kap * (zz[0] + zz[3])
    c = [s * s, kap * s / h**2, kap * s / (2 * h**4 * yy), kap * s / (6 * h**5 * yy**2), kap**2 / (24 * h**7 * yy**2)]
    X1, X2, X3 = ops
    X31 = mm(X3, X1)
    M = [[c[0] * int(i == j) + c[1] * X1[i][j] + c[2] * X2[i][j] + c[3] * X3[i][j] + c[4] * X31[i][j] for j in range(16)] for i in range(16)]
    return M
def rowvec_mat(r, M): return matvec(r, M)
I4 = eye()
fams = [
 ('W_a: exp(E14/3)', madd(I4, msc(Fr(1, 3), E(1, 4)))),
 ('torus x W_a: diag(2,1,1,3)(1+E14/5)', mm([[Fr(2), 0, 0, 0], [0, Fr(1), 0, 0], [0, 0, Fr(1), 0], [0, 0, 0, Fr(3)]], madd(I4, msc(Fr(1, 5), E(1, 4))))),
 ('R5: 1+2E24+3E34 (v=e1)', madd(I4, msc(2, E(2, 4)), msc(3, E(3, 4)))),
 ('R2: (1+2E12+3E13)(1+3E24-2E34) [p.c=0]', mm(madd(I4, msc(2, E(1, 2)), msc(3, E(1, 3))), madd(I4, msc(3, E(2, 4)), msc(-2, E(3, 4))))),
 ('R2 broken: (1+2E12+3E13)(1+E24+E34) [p.c=5]', mm(madd(I4, msc(2, E(1, 2)), msc(3, E(1, 3))), madd(I4, E(2, 4), E(3, 4)))),
 ('generic Borel element', g),
 ('W_a x R2: exp(E14/3)(1+2E12+3E13)(1+3E24-2E34)', mm(madd(I4, msc(Fr(1, 3), E(1, 4))), mm(madd(I4, msc(2, E(1, 2)), msc(3, E(1, 3))), madd(I4, msc(3, E(2, 4)), msc(-2, E(3, 4)))))),
]
vac = kron(vec(1, 0, 0, 0), vec(0, 0, 0, 1))
kap = Fr(3, 10)
Kl = Kop(kap, Lam, zz)
P('(6b) family | v=e1 g | u=e4 g^-T | closed form 5c | exact covariance <vac|rho(g)K[Lam] == <vac|K[G\']rho(g) (G\'=g Lam g^-1) | <vac|K_kappa[G\'] ~ <vac|rho(1-kappa G\') | (K_kappa[Lam] on v(x)u) ~ T_kappa | [G\'^n]_14 n=1,2,3 | off-diagonal part of G\'')
for name, gg in fams:
    gi_ = inv(gg); Gt = mm(mm(gg, Lam), gi_)        # Gt = g Lam g^-1
    Gp = tr(Gt)                                       # G' = g^-T Lam g^T = Gt^T   (found in (6a))
    v = [gg[0][j] for j in range(4)]
    u = [tr(gi_)[3][j] for j in range(4)]
    Rg_ = rhoM(gg)
    Kg = Kop(kap, Gp, zz)
    lhs = matvec(vac, Kg)
    a1 = matvec(kron(v, u), Kl) == matvec(matvec(vac, Kg), Rg_)
    one_m = madd(I4, msc(-kap, Gt)); tgt = kron([one_m[0][j] for j in range(4)], [tr(inv(one_m))[3][j] for j in range(4)])
    a2 = any(x != 0 for x in lhs) and rank([lhs, tgt]) == 1
    a3 = rank([matvec(kron(v, u), Kl), target(v, u, kap, zz)]) == 1
    Gn = [Gp, mm(Gp, Gp), mm(mm(Gp, Gp), Gp)]
    off = [[str(Gp[i][j]) if i != j and Gp[i][j] != 0 else '.' for j in range(4)] for i in range(4)]
    P(name, '|', [str(x) for x in v], [str(x) for x in u], '|', pred(v, u), '|', a1, '|', a2, '|', a3, '|', [str(G_[3][0]) for G_ in Gn], '|', ' '.join(''.join(f'{x:>6}' for x in r) for r in off))

# (c) Lie algebra h_r of pair-1 rotations in the span
def rhoE(i, j):
    # E_ij (x) 1 - 1 (x) E_ji acting on bras from the right
    A = E(i, j); B = E(j, i)
    return [[A[a][c] * int(b == d) - int(a == c) * B[b][d] for c in range(4) for d in range(4)] for a in range(4) for b in range(4)]
RHO = {(i, j): rhoE(i, j) for i in range(1, 5) for j in range(1, 5)}
def lie(v, u):
    r = kron(v, u); ops = pair_ops_z(zz); pw = words(v, u, ops)
    tg = [matvec(r, RHO[(i, j)]) for i in range(1, 5) for j in range(1, 5)]
    stab = nullspace(tg)
    ns = nullspace(tg + [[-x for x in w] for w in pw])
    hbasis, _ = rref([x[:16] for x in ns]); hbasis = [b for b in hbasis if any(x != 0 for x in b)]
    cur = list(stab); quo = []
    for x in hbasis:
        if rank(cur + [x]) > rank(cur): quo.append(x); cur.append(x)
    def xs(x):
        terms = []
        for i in range(4):
            for j in range(4):
                c = x[i * 4 + j]
                if c != 0: terms.append(f'{c}*E{i+1}{j+1}')
        return ' + '.join(terms)
    return rank(tg), len(stab), len(hbasis), len(quo), [xs(x) for x in quo], rank(pw)
def torus_test(v, u):
    ops = pair_ops_z(zz); pw = words(v, u, ops); rk = rank(pw)
    fs = {'Lambda': zz, 'Lambda^-1': [1 / z for z in zz], '1-kappa Lambda^-1 (3/10)': [1 - Fr(3, 10) / z for z in zz], 'Lambda^2': [z * z for z in zz], 'generic (2,5,5,-3)': [Fr(2), Fr(5), Fr(5), Fr(-3)]}
    res = {}
    for nm, f in fs.items():
        res[nm] = rank(pw + [kron([v[j] * f[j] for j in range(4)], [u[j] / f[j] for j in range(4)])]) == rk
    return res
refs = [('W_a (a=1/3)', vec(1, 0, 0, Fr(1, 3)), vec(Fr(-1, 3), 0, 0, 1)), ('R2: v=(1,2,3,0), u=(0,3,-2,5)', vec(1, 2, 3, 0), vec(0, 3, -2, 5)), ('R5: v=e1, u=(0,2,3,5)', vec(1, 0, 0, 0), vec(0, 2, 3, 5)),
        ('<0_w| limit (0,1,-1,0) (x) u_w', vec(0, 1, -1, 0), tw_u([Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)], zz)), ('R6: v=(0,2,3,0), u=(1,3,-2,4)', vec(0, 2, 3, 0), vec(1, 3, -2, 4)),
        ('random BPS (not admissible)', vec(1, 2, -1, 3), vec(Fr(5, 3), Fr(-1, 6), Fr(4, 3), 0))]
P('(9) reference | rank span | dim orbit tangent | dim stabiliser | dim h_r | dim h_r/stab | representatives of h_r mod stabiliser | torus-orbit test (f in span?)')
for nm, v, u in refs:
    rt, ds, dh, dq, reps, rs = lie(v, u)
    P(nm, '|', rs, '|', rt, '|', ds, '|', dh, '|', dq, '|', reps, '|', torus_test(v, u))
out.close()
