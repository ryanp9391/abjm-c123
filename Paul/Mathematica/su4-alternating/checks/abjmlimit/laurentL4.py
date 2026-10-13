import pickle
from taylor import *
from collections import Counter
p = P; F = ModP()
coef = pickle.load(open('dline_coef.pkl','rb'))['coef']
def is_reg(a,k): return a == 3 or (a == 2 and k <= 2)
def pair_terms(a,k): return [((a,k),1)] if not (a==2 and k<=2) else [((2,k),1),((1,k+1),-1)]
n = np.array([is_reg(w[0],w[1]) + is_reg(w[2],w[3]) for w in WIDX])
S = np.zeros((256,256), dtype=np.int64)
for j,w in enumerate(WIDX):
    for (q1,c1) in pair_terms(w[0],w[1]):
        for (q3,c3) in pair_terms(w[2],w[3]):
            S[j, WIDX.index((q1[0],q1[1],q3[0],q3[1]))] = (c1*c3) % p
SM = [ (S @ coef[m]) % p for m in range(len(coef)) ]
def Rm(m):
    return np.array([SM[m + n[j]][j] if m + n[j] < len(SM) else np.zeros(256, dtype=np.int64) for j in range(256)])
R0, R1, R2 = Rm(0), Rm(1), Rm(2)
print('rank R0', rank_mod(R0))
def inv_mod(A):
    n_ = A.shape[0]
    Rr, piv = rref_mod(np.concatenate([A % p, np.eye(n_, dtype=np.int64)], axis=1))
    return Rr[:, n_:]
R0i = inv_mod(R0)
assert not np.any((R0 @ R0i - np.eye(256, dtype=np.int64)) % p)
def vm(v, A): return (v @ A) % p
def bra(vs):
    out = np.array([1], dtype=np.int64)
    for v in vs: out = np.kron(out, np.array([F.c(x) for x in v], dtype=np.int64)) % p
    return out
tests = {'product bra (generic)': bra([[1,Fr(2,3),Fr(-1,2),Fr(3,5)],[Fr(2,7),-1,Fr(1,3),1],[Fr(-4,3),1,Fr(5,2),Fr(1,9)],[1,Fr(3,4),Fr(-2,5),Fr(7,3)]]),
         'BPS v w v w (v.w=0)': None}
v = [1, Fr(-2), 3, 1]; r = [2, 1, -1, 3]; vv = sum(Fr(x)*x for x in v); vr = sum(Fr(a)*b for a,b in zip(v,r))
w = [Fr(x) - vr/vv*Fr(y) for x,y in zip(r,v)]
tests['BPS v w v w (v.w=0)'] = bra([v,w,v,w])
for name, b in tests.items():
    b0 = vm(b, R0i); b1 = (-vm(vm(vm(b, R0i), R1), R0i)) % p
    b2 = (vm(vm(vm(vm(vm(b,R0i),R1),R0i),R1),R0i) - vm(vm(vm(b,R0i),R2),R0i)) % p
    beta = [b0, b1, b2]
    # Laurent coefficients of c: c = sum_j beta_j(eps) eps^{-n_j} S_j
    L = {}
    for order in (-2, -1, 0):
        acc = np.zeros(256, dtype=np.int64)
        for j in range(256):
            m = order + n[j]
            if 0 <= m <= 2 and beta[m][j]:
                acc = (acc + beta[m][j] * S[j]) % p
        L[order] = acc
    pole = []
    for x in range(256):
        po = 0
        for order in (-2,-1):
            if L[order][x]: po = -order; break
        pole.append(po)
    pc = Counter(pole)
    # which words have eps^-2: classify per pair
    def cls(a,k): return 'S' if (a==3 or (a==2 and k<=2) or (a==1 and 1<=k<=3)) else 'R'
    pattern = Counter((cls(w[0],w[1]), cls(w[2],w[3]), pole[i]) for i,w in enumerate(WIDX))
    print(name, ': pole-order counts', sorted(pc.items()), ' by per-pair class (S=singular set, R=regular):', sorted(pattern.items()))
    # finite combinations check: residues at eps^-2 equal beta0_j for double-reg rows
