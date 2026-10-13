# Same-point generating functions P(u0) = sum_k c_k t_k(u0 + s_k): which coefficient vectors c make it localise on
# pair 1 for several references simultaneously (generic theta and ABJM point; shifts s_k in {0, h, 2h, 3h, -h})
import numpy as np, itertools
from pgbase import *
h = Fr(1,2)
def refs_and_chain(abjm, zz):
    thF1, thF2 = Fr(-3,5), Fr(-11,8)
    th = [thF1 + (2*h if abjm else Fr(2,7)+Fr(3,5)), thF1, thF2 + (2*h if abjm else Fr(8,11)+Fr(11,8)), thF2]
    ch = Chain(['A','F','A','F'], th, h, diagm(zz))
    return ch, th
z = [Fr(3,2), Fr(5,4), Fr(2,3), Fr(4,7)]
rng = np.random.default_rng(9)
def Wp(a): return [np.array([1,0,0,a], float), np.array([-a,0,0,1], float)]
v = rng.normal(size=4); u = rng.normal(size=4); u -= (v@u)/(v@v)*v
refs = [Wp(1/3), Wp(-0.4), [v, u]]
def nonlocal_part(out, rest):
    M = out.reshape(16, 16); ps = M@rest/(rest@rest)
    return (M - np.outer(ps, rest)).reshape(-1)
for abjm in (False, True):
    ch, th = refs_and_chain(abjm, z)
    print('==== ABJM point' if abjm else '==== generic theta')
    for base, bname in ((th[1], 'thF'), (th[0], 'thA')):
        for m in range(-3, 5):
            u0 = base + m*h
            # all k at the same point (Talalaev-type) and the 'last column at u0' variant u0 + (k-1)h
            for variant in ('same', 'stair'):
                cols = []
                for k in range(0, 5):
                    uk = u0 if variant == 'same' else u0 + (k-1)*h
                    blocks = []
                    for (vv, uu) in refs:
                        ref = prod_state([vv, uu, vv, uu], dtype=float); rest = np.kron(vv, uu)
                        out = (ch.apply_t(k, uk, ref, reg=True) if k > 0 else ref).astype(float)
                        blocks.append(nonlocal_part(out, rest))
                    cols.append(np.concatenate(blocks))
                A = np.array(cols).T
                sv = np.linalg.svd(A, compute_uv=False)
                nul = 5 - np.sum(sv > 1e-9*max(sv[0], 1e-300))
                # individually localising k's
                indiv = [k for k in range(5) if np.linalg.norm(A[:, k]) < 1e-9*max(1, np.abs(A).max())]
                if nul > len(indiv):
                    print('  u0 = %s%+dh (%s): localising combinations dim %d > individually localising %s  <-- NEW' % (bname, m, variant, nul, indiv))
    print('  (scan done; only lines marked NEW would indicate a localising combination beyond single words)')
