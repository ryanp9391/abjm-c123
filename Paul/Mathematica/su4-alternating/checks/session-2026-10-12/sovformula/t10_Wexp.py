"""SoV expansion of HW-orbit references in the <0|-word basis (sov_l4 standard point, generic theta, exact mod p):
<W_a^L| (diag frame (e1 + a e4)(x)(e4 - a e1) per pair) mapped to the companion frame; count nonzero coefficients."""
import pickle, numpy as np, sys
from fractions import Fraction as Fr
sys.path.insert(0, '../sovq')
import sov_l4 as S
from core import modinv
p = 2147483629 // 32 * 0 + 67108859
W = pickle.load(open('../sovq/cache_sov_l4/words.pkl', 'rb'))     # 256 x 256 Fractions, rows = words
def modarr(M):
    return np.array([[(x.numerator % p) * pow(x.denominator % p, -1, p) % p for x in row] for row in M], dtype=np.int64)
Wm = modarr(W)
Wi = modinv(Wm, p)
z = S.Z0
V = S.vandermonde(z); Vi = S.minv(V); VT = S.mT(V)
def comp_pair(vd, ud):
    vc = [sum(vd[i] * VT[i][j] for i in range(4)) for j in range(4)]          # v_c = v_d V^T
    uc = [sum(ud[i] * Vi[i][j] for i in range(4)) for j in range(4)]          # u_c = u_d V^{-1}
    return vc, uc
def ratrec(a):
    from core import ratrec as rr
    return rr(a, p)
for a in (Fr(0), Fr(-1), Fr(1, 3)):
    vd = [Fr(1), 0, 0, a]; ud = [-a, 0, 0, Fr(1)]
    vc, uc = comp_pair(vd, ud)
    bra = S.product_vec([vc, uc, vc, uc])
    b = modarr([bra])[0]
    c = (b @ Wi) % p if False else None
    # c W = bra  ->  c = bra W^{-1}
    c = np.zeros(256, dtype=np.int64)
    for j in range(256):
        c = (c + int(b[j]) * Wi[j]) % p if b[j] else c
    nz = [(S.WORDS[i], c[i]) for i in range(256) if c[i] != 0]
    print('a =', a, ': nonzero SoV coefficients of <W_a^L| :', len(nz))
    if len(nz) <= 30:
        for wd, x in nz: print('    ', wd, ratrec(x))
    # by-word statistics: which labels occur
    labs = [set(wd[s] for wd, _ in nz) for s in range(4)]
    print('   labels used per site:', labs)
