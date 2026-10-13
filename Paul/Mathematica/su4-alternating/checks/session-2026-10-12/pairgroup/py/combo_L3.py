# L = 3 pairs: do linear combinations of words at pair-1 points (singles t_k(theta_2+m h), products of two) give
# W_{a'} (x) W_a (x) W_a universally (3 values of a)?  ABJM point; z2 != z3 or z2 = z3.
import numpy as np, itertools, sys, time
from pgbase import *
h = Fr(1,2)
zz = [Fr(3,2), Fr(5,4), Fr(5,4), Fr(4,7)] if sys.argv[1] == 'z23' else [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
abjm = sys.argv[2] == 'abjm' if len(sys.argv) > 2 else True
thF = [Fr(-3,5), Fr(-11,8), Fr(-7,9)]
th = []
for i, f in enumerate(thF): th += [f + (2*h if abjm else [Fr(31,35), Fr(39,40), Fr(13,9)][i]), f]
ch = Chain(['A','F']*3, th, h, diagm(zz))
def Wp(a): return [np.array([1,0,0,a], float), np.array([-a,0,0,1], float)]
avals = (1/3, -0.4, 1.7)
singles = [(k, m) for k in (1,2,3) for m in range(-2, 5)]
t0 = time.time()
Fcols = {(): [], }; NLcols = {(): []}
outs = {}
for r, a in enumerate(avals):
    ref = prod_state(Wp(a)*3, dtype=float); rest = prod_state(Wp(a)*2, dtype=float).reshape(-1)
    outs[(r, ())] = ref
    for s in singles:
        outs[(r, (s,))] = ch.apply_t(s[0], th[1] + s[1]*h, ref, reg=True).astype(float)
    for s1, s2 in itertools.combinations_with_replacement(singles, 2):
        outs[(r, (s1, s2))] = ch.apply_t(s2[0], th[1] + s2[1]*h, outs[(r, (s1,))], reg=True).astype(float)
print('applications done %.0f s' % (time.time()-t0))
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
scale = np.maximum(np.abs(F).max(0), np.abs(NL).max(0)); scale[scale == 0] = 1
F = F/scale; NL = NL/scale
u_, s_, vh = np.linalg.svd(NL, full_matrices=True)
rk = np.sum(s_ > 1e-9*s_[0]); null = vh[rk:].T
img = F @ null
si = np.linalg.svd(img, compute_uv=False); rank_img = np.sum(si > 1e-8*si[0])
z = [float(x) for x in zz]; kap = 0.3
tg = np.concatenate([(1-kap*z[0])**2*prod_state(Wp(a*(1-kap*z[3])/(1-kap*z[0])), dtype=float).reshape(-1) for a in avals])
c, *_ = np.linalg.lstsq(img, tg, rcond=None)
print('L=3', 'ABJM' if abjm else 'generic', 'z', [str(x) for x in zz], ': #words', len(words), ' localising-combination dim', null.shape[1],
      ' rank of localised pair states', rank_img, ' target residual %.2e' % (np.linalg.norm(img@c - tg)/np.linalg.norm(tg)))
# sparse localising combinations: pairs / triples of words whose non-local parts are linearly dependent
NLn = NL/np.maximum(np.linalg.norm(NL, axis=0), 1e-300)
nz = np.linalg.norm(NL, axis=0) > 1e-9
idx = np.where(nz)[0]
Gm = NLn[:, idx].T @ NLn[:, idx]
found = []
for p in range(len(idx)):
    for q in range(p+1, len(idx)):
        if abs(abs(Gm[p, q]) - 1) < 1e-10:
            found.append((idx[p], idx[q], Gm[p, q]))
print('pairs of words with proportional non-local parts:', len(found))
# localised outputs of these 2-word combos; which are new w.r.t. span of individually localising words
indiv = [i for i in range(len(words)) if not nz[i]]
B0 = F[:, indiv]
def newdir(vec):
    c, *_ = np.linalg.lstsq(B0, vec, rcond=None); return np.linalg.norm(B0@c - vec)/max(np.linalg.norm(vec), 1e-300)
for (i, j, g) in found[:400]:
    lam_ = -(np.linalg.norm(NL[:, i])/np.linalg.norm(NL[:, j]))*np.sign(g)
    vec = F[:, i] + lam_*F[:, j]
    if np.linalg.norm(vec) > 1e-8 and newdir(vec) > 1e-6:
        print('  NEW localised direction from', words[i], '+ c *', words[j], ' (rel. distance from individual span %.2e)' % newdir(vec))
