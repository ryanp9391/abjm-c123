# Do LINEAR COMBINATIONS of (non-individually-localising) words localise on pair 1, universally (3 references at once)?
# columns: 1, single t_k(theta_2 + m h) (k=1,2,3; m=-2..4, reg), and all products of two.  4 sites, ABJM point or generic.
import numpy as np, itertools, sys
from pgbase import *
h = Fr(1,2)
abjm = (sys.argv[1] == 'abjm') if len(sys.argv) > 1 else True
zz = [Fr(3,2), Fr(5,4), Fr(5,4), Fr(4,7)] if (len(sys.argv) > 2 and sys.argv[2] == 'z23') else [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
thF1, thF2 = Fr(-3,5), Fr(-11,8)
th = [thF1 + (2*h if abjm else Fr(31,35)), thF1, thF2 + (2*h if abjm else Fr(11,8)+Fr(8,11)), thF2]
ch = Chain(['A','F','A','F'], th, h, diagm(zz))
rng = np.random.default_rng(9)
def Wp(a): return [np.array([1,0,0,a], float), np.array([-a,0,0,1], float)]
refsel = sys.argv[3] if len(sys.argv) > 3 else 'mixed'
if refsel == 'W':
    refs = [Wp(1/3), Wp(-0.4), Wp(1.7)]
else:
    v = rng.normal(size=4); u = rng.normal(size=4); u -= (v@u)/(v@v)*v
    refs = [Wp(1/3), Wp(-0.4), [v, u]]
singles = [(k, m) for k in (1,2,3) for m in range(-2, 5)]
words = [()] + [(s,) for s in singles] + list(itertools.combinations_with_replacement(singles, 2))
F = []; NL = []
cache = {}
for w in words:
    fparts = []; nparts = []
    for r, (vv, uu) in enumerate(refs):
        ref = prod_state([vv, uu, vv, uu], dtype=float); rest = np.kron(vv, uu)
        psi = ref
        for (k, m) in w:
            key = (r,) + tuple(w[:w.index((k,m))+1])
            psi = ch.apply_t(k, th[1] + m*h, psi, reg=True).astype(float)
        M = psi.reshape(16, 16); ps = M@rest/(rest@rest)
        fparts.append(psi.reshape(-1)); nparts.append((M - np.outer(ps, rest)).reshape(-1))
    F.append(np.concatenate(fparts)); NL.append(np.concatenate(nparts))
F = np.array(F).T; NL = np.array(NL).T
scale = np.abs(F).max(axis=0); scale[scale == 0] = 1
F = F/scale; NL = NL/scale
u_, s_, vh = np.linalg.svd(NL, full_matrices=True)
r = np.sum(s_ > 1e-9*s_[0])
null = vh[r:].T
img = F @ null
si = np.linalg.svd(img, compute_uv=False)
rank_img = np.sum(si > 1e-8*si[0]) if len(si) else 0
indiv = [i for i in range(len(words)) if np.linalg.norm(NL[:, i]) < 1e-9]
si2 = np.linalg.svd(F[:, indiv], compute_uv=False)
rank_ind = np.sum(si2 > 1e-8*si2[0])
print('ABJM' if abjm else 'generic', 'z =', [str(x) for x in zz], 'refs:', refsel)
print('  #words', len(words), ' localising-combination space dim', null.shape[1], ' rank of localised outputs', rank_img,
      ' | individually localising words', len(indiv), 'rank', rank_ind)
# target test: one combination (same coefficients for all refs) giving (1-kz1)^2 W_{a'} (x) W_a on each W ref
if refsel == 'W':
    z = [float(x) for x in zz]; kap = 0.3
    tg = []
    for (vv, uu), a in zip(refs, (1/3, -0.4, 1.7)):
        ap = a*(1-kap*z[3])/(1-kap*z[0])
        tg.append((1-kap*z[0])**2*prod_state(Wp(ap) + Wp(a), dtype=float).reshape(-1))
    tg = np.concatenate(tg)
    c, *_ = np.linalg.lstsq(img, tg, rcond=None)
    print('  target W_{a\'} (x) W_a for all three a with one combination: rel. residual %.2e' % (np.linalg.norm(img@c - tg)/np.linalg.norm(tg)))
