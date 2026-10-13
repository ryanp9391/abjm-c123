"""2603-type 'twist-dressed' references at generic theta (sov_l4 standard point, L = 2 pairs, exact mod p):
<W| = <0| prod_A rho_A-site(gA(Gc^T)) prod_F (u -> u gF(Gc)), i.e. diag frame A: lambda_i gA(z_i), F: gF(z_j).
Count nonzero SoV coefficients in the <0|-word basis."""
import pickle, numpy as np, sys, itertools
from fractions import Fraction as Fr
import sov_l4 as S
from core import modinv
p = 67108859
W = pickle.load(open('../sovq/cache_sov_l4/words.pkl', 'rb'))
def m(x): x = Fr(x); return (x.numerator % p) * pow(x.denominator % p, -1, p) % p
Wi = modinv(np.array([[m(x) for x in r] for r in W], dtype=np.int64), p)
z = S.Z0; Gc = S.Gc(z); GcT = S.mT(Gc)
I4 = [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
def inv(A): return S.minv(A)
def lin(a, A, b, B): return [[a * A[i][j] + b * B[i][j] for j in range(4)] for i in range(4)]
def vm(v, M): return [sum(v[i] * M[i][j] for i in range(4)) for j in range(4)]
e1 = [Fr(1), 0, 0, 0]; e4 = [0, 0, 0, Fr(1)]
def expand(vA1, uF2, vA3, uF4):
    b = S.product_vec([vA1, uF2, vA3, uF4])
    bm = np.array([m(x) for x in b], dtype=np.int64)
    c = np.zeros(256, dtype=np.int64)
    for j in np.nonzero(bm)[0]: c = (c + int(bm[j]) * Wi[j]) % p
    return [(S.WORDS[i], int(c[i])) for i in range(256) if c[i] != 0]
for om in (Fr(2, 7), Fr(-3, 4)):
    RA = inv(lin(1, I4, -om, GcT)); RF = inv(lin(1, I4, -om, Gc))
    cases = {'A sites (1-w Gc^T)^-1': (vm(e1, RA), e4, vm(e1, RA), e4),
             'F sites (1-w Gc)^-1': (e1, vm(e4, RF), e1, vm(e4, RF)),
             'A site 1 only': (vm(e1, RA), e4, e1, e4), 'F site 2 only': (e1, vm(e4, RF), e1, e4),
             'A site 3 only': (e1, e4, vm(e1, RA), e4), 'F site 4 only': (e1, e4, e1, vm(e4, RF)),
             'all sites': (vm(e1, RA), vm(e4, RF), vm(e1, RA), vm(e4, RF)),
             'A: (1-wGc^T), F: (1-wGc)^-1 ...': (vm(e1, lin(1, I4, -om, GcT)), vm(e4, RF), vm(e1, lin(1, I4, -om, GcT)), vm(e4, RF))}
    for name, vs in cases.items():
        nz = expand(*vs)
        labs = [sorted(set(w[s] for w, _ in nz)) for s in range(4)]
        print(f'omega={om}  {name:35s} nonzero words: {len(nz):3d}   labels per site: {labs}')
