# Task 1f: on the 2-pair block (L = 2 chain = reduced chain of any L by block localisation):
# which W_b (x) W_c lie in the 25-dim span of multilinear monomials?  Intersection of span{u_i (x) u_j} with the monomial span.
import numpy as np, itertools
from fractions import Fraction as Fr
import eng, lp, modlin
from eng import P, modp
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
thF = [Fr(-3, 5), Fr(-11, 8)]
ch = lp.chain_abjm(2, thF, h, z, 'p')
for a in (Fr(1, 3), Fr(-1)):
    ref = lp.Wstate(ch, [a, a])
    vecs = {(): ref}
    for p in (1, 2):
        new = {}
        for key, v in vecs.items():
            new[key + (0,)] = v
            for i, x in enumerate(lp.Xs(ch, v, p)): new[key + (i + 1,)] = x
        vecs = new
    mons = [vecs[k] for k in sorted(vecs)]
    e = np.eye(4, dtype=np.int64)
    u = [np.multiply.outer(e[0], e[3]), (np.multiply.outer(e[3], e[3]) - np.multiply.outer(e[0], e[0])) % P, (-np.multiply.outer(e[3], e[0])) % P]
    U = [np.multiply.outer(u[i], u[j]).reshape(-1) % P for i in range(3) for j in range(3)]
    # intersection: solve sum_m c_m mon_m - sum_ij d_ij U_ij = 0 -> nullspace
    cols = [m.reshape(-1) for m in mons] + [(-x) % P for x in U]
    rk, cons, sol, null = modlin.solve(cols, np.zeros_like(cols[0]))
    print('a = %s: rank of [monomials | 9 W-tensors] = %d of %d; intersection dim = %d' % (a, rk, len(cols), len(null)))
    for nv in null:
        d = [modlin.ratrec(int(x)) for x in nv[25:]]
        print('    d_ij (coeffs of u_i u_j, i,j=0..2):', [str(x) for x in d])
