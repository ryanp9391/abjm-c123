import time, pickle
from taylor import *
from itertools import combinations
F = ModP(); h = H; p = F.p
tF1, tF2 = Fr(-3,5), Fr(-11,8)
th = [tF1+2*h, tF1, tF2+2*h, tF2]
ch = Chain(F, TYPES, th, h, G)
def Bop(u):
    u = Fr(u)
    Tu = ch.T(u)[0]
    acc = F.zeros(256)
    for k1 in (1,2,3):
        for K2 in combinations((1,2,3), 2):
            A = Tu[k1-1][0]
            Bm = ch.qmin(tuple(k-1 for k in K2), (0, k1), u+h, twisted=False)
            C = ch.qmin((0,1,2), (0, K2[0], K2[1]), u+2*h, twisted=False)
            acc = (acc + F.mm(F.mm(A, Bm), C)) % p
    return acc
t0 = time.time()
B1 = Bop(Fr(37,100)); 
print('B built', time.time()-t0, flush=True)
R = pickle.load(open('regL4.pkl','rb'))['R']
def is_reg(a,k): return a == 3 or (a == 2 and k <= 2)
# eigen test for each row: rank{v, vB} == 1 ?
eig = []
for i,w in enumerate(WIDX):
    v = R[i]; vB = v @ B1 % p
    eig.append(rank_mod(np.array([v, vB])) == 1)
eig = np.array(eig)
regtype = np.array([(is_reg(w[0],w[1]) and w[0]==2) or (is_reg(w[2],w[3]) and w[2]==2) for w in WIDX])
print('# regularised rows that are B-eigen:', eig.sum(), 'of 256; non-eigen rows all contain a W (a=2 regularised) factor:', set(regtype[~eig]), ' count', (~eig).sum())
# Jordan: for non-eigen rows v, find beta s.t. vB - beta v lies in span of eigen rows with eigenvalue beta
# compute eigenvalue of eigen rows
def eigval(v, vB):
    j = np.nonzero(v)[0][0]; return int(vB[j]) * pow(int(v[j]), p-2, p) % p
betas = {}
for i in np.nonzero(eig)[0]:
    v = R[i]; betas[i] = eigval(v, v @ B1 % p)
from collections import Counter
cnt = Counter(betas.values())
print('# distinct B eigenvalues among eigen rows:', len(cnt), ' multiplicity pattern', sorted(Counter(cnt.values()).items()))
# Jordan test on one non-eigen row
jcount = 0
for i in np.nonzero(~eig)[0]:
    v = R[i]; vB = v @ B1 % p
    ok = False
    for b in set(betas.values()):
        r = (vB - b*v) % p
        same = [R[j] for j in betas if betas[j]==b]
        if rank_mod(np.array(same + [r])) == len(same):
            ok = True; break
    jcount += ok
print('# non-eigen regularised rows with vB = beta v + (eigen rows of same beta):', jcount)
# rank of (B - beta)^2 vs (B - beta) for one degenerate beta
# length-3 chains: rows with two W factors
both = [i for i,w in enumerate(WIDX) if (w[0]==2 and w[1]<=2) and (w[2]==2 and w[3]<=2)]
ok3 = 0
groups = {}
for i in np.nonzero(eig)[0]: groups.setdefault(betas[i], []).append(i)
# assign beta to single-W rows
for i in np.nonzero(~eig)[0]:
    if i in both: continue
    v = R[i]; vB = v @ B1 % p
    for b, idx in groups.items():
        r = (vB - b*v) % p
        if rank_mod(np.array([R[j] for j in idx] + [r])) == len(idx):
            betas[i] = b; break
for i in both:
    v = R[i]; vB = v @ B1 % p
    for b in set(betas.values()):
        r = (vB - b*v) % p
        same = [R[j] for j in betas if betas[j]==b]
        if rank_mod(np.array(same + [r])) == rank_mod(np.array(same)):
            ok3 += 1; break
print('# double-W rows with vB = beta v + (same-beta eigen & single-W rows):', ok3, 'of', len(both))
# Jordan block sizes: for each beta, dim of generalized eigenspace (# rows assigned) 
sizes = Counter(Counter(betas.values()).values())
print('generalised-eigenspace dims (single-W rows assigned):', sorted(sizes.items()))
cls = {}
for b, c in Counter(betas.values()).items(): cls.setdefault(c, []).append(b)
I = np.eye(256, dtype=np.int64)
for c in sorted(cls):
    b = cls[c][0]
    if c == 3: c_true = 4
    else: c_true = c
    Bb = (B1 - b*I) % p
    r1 = rank_mod(Bb); Bb2 = F.mm(Bb,Bb); r2 = rank_mod(Bb2); r3 = rank_mod(F.mm(Bb2,Bb))
    print('class with', c_true, 'basis rows: dim ker(B-b)=', 256-r1, ' dim ker(B-b)^2=', 256-r2, ' dim ker (B-b)^3=', 256-r3)
