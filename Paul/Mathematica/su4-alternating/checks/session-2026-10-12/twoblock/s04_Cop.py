"""Search for a GLMR-type C-operator (commuting family, untwisted T) and test separability of <B|y> on its right eigenvectors."""
import numpy as np, itertools, time
from core import *
from itertools import permutations, combinations
R = Ring('float'); I4 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
ch = Chain(R, G=I4)    # untwisted
h = float(h0)
def T(u):
    return ch.mono(Fr(u).limit_denominator(10**9) if not isinstance(u, Fr) else u)
Tc = {}
def Tm(u):
    key = u
    if key not in Tc: Tc[key] = ch.mono(u)
    return Tc[key]
def qmin(rows, cols, v, sgn_shift=-1, rev=False):
    """sum_sigma sgn T_{rows[s(1)],cols[1]}(v) T_{rows[s(2)],cols[2]}(v + sgn_shift h) ... ; rev: reverse product order"""
    a = len(rows); out = 0
    for sg in permutations(range(a)):
        s = perm_sign(sg)
        facs = [Tm(v + sgn_shift * k * h0)[rows[sg[k]]][cols[k]] for k in range(a)]
        if rev: facs = facs[::-1]
        P = facs[0]
        for f in facs[1:]: P = P @ f
        out = out + s * P
    return out
def B_RV(u):
    out = 0
    for k1 in range(3):
        for K2 in combinations(range(3), 2):
            out = out + Tm(u)[k1][0] @ qmin(list(K2), [0, k1 + 1], u + h0) @ qmin([0, 1, 2], [0, K2[0] + 1, K2[1] + 1], u + 2 * h0)
    return out
u1, u2 = Fr(37, 100), Fr(-113, 100)
Bu, Bv = B_RV(u1), B_RV(u2)
print('RV B commutes:', np.abs(Bu @ Bv - Bv @ Bu).max() / np.abs(Bu).max() ** 2)
# candidate C operators: reflection i -> 3-i (0-based), optional transpose, shift sign, product order
def refl(i): return 3 - i
def C_cand(u, transpose, sh, order, rev):
    out = 0
    for k1 in range(3):
        for K2 in combinations(range(3), 2):
            r1 = [refl(k1)]; c1 = [refl(0)]
            r2 = [refl(k) for k in K2]; c2 = [refl(0), refl(k1 + 1)]
            r3 = [refl(0), refl(1), refl(2)]; c3 = [refl(0), refl(K2[0] + 1), refl(K2[1] + 1)]
            if transpose: (r1, c1), (r2, c2), (r3, c3) = (c1, r1), (c2, r2), (c3, r3)
            if transpose and len(r1) != len(c1): pass
            f1 = qmin(r1, c1, u, -1 if sh < 0 else 1, rev) if len(r1) == len(c1) else None
            f2 = qmin(r2, c2, u + sh * h0, -1 if sh < 0 else 1, rev)
            f3 = qmin(r3, c3, u + 2 * sh * h0, -1 if sh < 0 else 1, rev)
            fs = [f1, f2, f3] if order == 0 else [f3, f2, f1]
            out = out + fs[0] @ fs[1] @ fs[2]
    return out
results = []
for transpose in (False, True):
    for sh in (1, -1):
        for order in (0, 1):
            for rev in (False, True):
                try:
                    Cu, Cv = C_cand(u1, transpose, sh, order, rev), C_cand(u2, transpose, sh, order, rev)
                except Exception as ex:
                    print('fail', ex); continue
                cm = np.abs(Cu @ Cv - Cv @ Cu).max() / max(np.abs(Cu).max() ** 2, 1e-300)
                ev = np.linalg.eigvals(Cu)
                d = np.abs(ev[:, None] - ev[None, :]) + np.eye(256) * 1e300
                print(f'transpose={transpose} shift={sh:+d} order={order} rev={rev}: rel commutator {cm:.2e}, |C|max {np.abs(Cu).max():.2e}, min eig gap {d.min():.2e}, rank {np.linalg.matrix_rank(Cu)}')
