"""Candidate test: product over pairs of one-pair rotated overlaps (pair 1 rotated / pair 2 at SoV vacuum, and vice versa), mp 50 digits on all 256 states.
Also export the exact SoV coefficients c^(n)_x of <W_a|R^-1 and <0_w| to JSON."""
from wa_sov_core import *
import numpy as np, json
import sov_l4 as S
import mpmath as mp
mp.mp.dps = 50
eig = cached('eigen', S.eigenstates)
sts = [S.Eig(d) for d in eig]
right, left = cached('lr', lambda: S.companion_vectors(eig))
cn = pickle.load(open(os.path.join(CACHE, 'cn.pkl'), 'rb'))
V = vandermonde(Z0); Vi = minv(V)
vz = [[V[i][j] for i in range(4)] for j in range(4)]; ut = [Vi[j] for j in range(4)]
a0 = Fr(3, 7)
Av = [vz[0][i] + a0 * vz[3][i] for i in range(4)]; Fu = [ut[3][i] - a0 * ut[0][i] for i in range(4)]
e1, e4 = unit(0), unit(3)
full = product_vec([Av, Fu, Av, Fu]); p1 = product_vec([Av, Fu, e1, e4]); p2 = product_vec([e1, e4, Av, Fu]); vac = product_vec([e1, e4, e1, e4])
def ov(b, k): return mp.fsum(S.tomp(b[c]) * right[k][c] for c in range(256) if b[c] != 0)
rows = []
for k, st in enumerate(sts):
    f = ov(full, k)
    if abs(f) < mp.mpf(10) ** -30: continue
    o0 = right[k][51]
    rows.append((st.weight, f / o0, ov(p1, k) / o0, ov(p2, k) / o0))
print('states with nonzero <W_a|Psi>:', len(rows))
ratio = [r[1] / (r[2] * r[3]) for r in rows]
print('G / (P1 P2): min |.| %s max |.| %s  (a constant ratio would signal a product structure)' % (mp.nstr(min(abs(x) for x in ratio), 6), mp.nstr(max(abs(x) for x in ratio), 6)))
# the same for the sum form: is G = c1*P1 + c2*P2 + c0 (single-pair rotations) ?
import itertools
A = mp.matrix([[r[2], r[3], 1] for r in rows]); bvec = mp.matrix([r[1] for r in rows])
sol = mp.lu_solve(A, bvec)   # least squares
res = A * sol - bvec
print('least-squares G = c1 P1 + c2 P2 + c0: residual norm / |G| =', mp.nstr(mp.norm(res) / mp.norm(bvec), 5))
json.dump({'meaning': 'c^(n)_x = coefficient of a^n of the SoV expansion <W_a|R^-1 = sum_x c_x(a) <x|, x = (a1,k2,a3,k4) in itertools.product order; entries "num/den"',
           'words': [list(w) for w in WORDS], 'components': {str(n): [str(c) for c in cn[n]] for n in range(5)}},
          open(os.path.join(CACHE, '..', 'Wa_sov_coefficients.json'), 'w'))
print('exported')
