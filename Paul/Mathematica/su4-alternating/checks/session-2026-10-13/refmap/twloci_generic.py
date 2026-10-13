"""twloci_generic.py -- twisted vacuum at GENERIC z (z2 != z3): which loci in w admit a kappa-family?  (exact)"""
from helpers import *
import itertools
zg = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
wgen = [Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)]
def run(name, w):
    v, u = tw_form(w, zg)
    rk, rkT, ok, cs = test(v, u, zg)
    pis = [v[j] * u[j] for j in range(4)]
    print(name, '| v =', [str(x) for x in v], '| u =', [str(x) for x in u], '| pi_j =', [str(x) for x in pis], '| rank', rk, rkT, '| K-family:', ok)
run('generic w', wgen)
for a in range(4):
    w = list(wgen); w[0] = zg[a]; run('w1 = z%d' % (a + 1), w)
for a, b in itertools.combinations(range(4), 2):
    w = list(wgen); w[1] = zg[a]; w[2] = zg[b]; run('w2 = z%d, w3 = z%d' % (a + 1, b + 1), w)
w = list(wgen); w[1] = zg[1]; w[2] = zg[2]; w[0] = Fr(7, 5); run('w2,w3 = z2,z3; w1 generic', w)
w = [Fr(7, 5), zg[0], zg[1], Fr(5, 9)]; run('w2,w3 = z1,z2; w1 generic', w)
w = list(wgen); w[3] = zg[3]; run('w4 = z4 (u independent of w4)', w)
