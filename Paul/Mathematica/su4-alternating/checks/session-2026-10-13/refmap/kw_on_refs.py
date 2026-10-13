"""kw_on_refs.py -- does the fixed-coefficient W_a operator K_kappa act as the common rotation on the other admissible classes?  And: R6 closed-form operator.  (python fallback)"""
from helpers import *
y = Fr(5, 4); zz = [Fr(3, 2), y, y, Fr(4, 7)]
ops = pair_ops_z(zz)
def Kcoef_W(kap):
    s = 1 - kap * (zz[0] + zz[3])
    return [s * s, kap * s / h**2, kap * s / (2 * h**4 * y), kap * s / (6 * h**5 * y**2), kap**2 / (24 * h**7 * y**2)]
def Kcoef_R6(kap):
    e1 = zz[0] + y + zz[3]; e2 = y * zz[0] + y * zz[3] + zz[0] * zz[3]
    return [1 - kap * e1 + kap**2 * e2, 2 * kap, 8 * kap**2, Fr(0), Fr(0)]
refs = [('W_a', vec(1, 0, 0, Fr(1, 3)), vec(Fr(-1, 3), 0, 0, 1)),
        ('<0|comp limit (0,1,-1,0)(x)(1,1,1,1)', vec(0, 1, -1, 0), vec(1, 1, 1, 1)),
        ('<0_w| limit (0,1,-1,0)(x)u_w', vec(0, 1, -1, 0), tw_u([Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)], zz)),
        ('R6 (0,2,3,0)(x)(1,3,-2,4)', vec(0, 2, 3, 0), vec(1, 3, -2, 4)),
        ('R5 e1 (x) (0,2,3,5)', vec(1, 0, 0, 0), vec(0, 2, 3, 5)),
        ('R2 (1,2,3,0)(x)(0,3,-2,5)', vec(1, 2, 3, 0), vec(0, 3, -2, 5))]
for kap in (Fr(3, 10), Fr(-7, 5), Fr(7, 4)):
    for nm, v, u in refs:
        pw = words(v, u, ops); T = [t / (1 - kap * y) ** 2 for t in target(v, u, kap, zz)]
        res = []
        for lab, cf in (('K_W', Kcoef_W(kap)), ('K_R6', Kcoef_R6(kap))):
            lhs = [sum(c * w[j] for c, w in zip(cf, pw)) for j in range(16)]
            res.append((lab, lhs == T))
        print('kappa', kap, nm, res)
