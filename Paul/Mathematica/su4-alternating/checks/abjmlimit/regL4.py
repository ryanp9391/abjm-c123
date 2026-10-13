import time, pickle
from taylor import *
from collections import Counter
F = ModP(); h = H; p = F.p
tF1, tF2 = Fr(-3,5), Fr(-11,8)
th = [tF1+2*h, tF1, tF2+2*h, tF2]
ch = Chain(F, TYPES, th, h, G)
d = pickle.load(open('taylor_pbp.pkl','rb')); TA1, TA3 = d['TA1'], d['TA3']
Y2 = ch.t(1, th[1]); Y4 = ch.t(1, th[3])
vac = vacuum(TYPES, F)
W1 = (TA1[2][1] - F.mm(TA1[1][1], Y2)) % p; Z1 = TA1[3][1]
W3 = (TA3[2][1] - F.mm(TA3[1][1], Y4)) % p; Z3 = TA3[3][1]
def op(TA, W, Z, a, k):
    if a < 2: return TA[a][0]
    if a == 2: return W if k <= 2 else TA[2][0]
    return Z
def regwords():
    rows = []
    for (a1,k2,a3,k4) in WIDX:
        v = vac @ op(TA1,W1,Z1,a1,k2) % p
        for _ in range(k2): v = v @ Y2 % p
        v = v @ op(TA3,W3,Z3,a3,k4) % p
        for _ in range(k4): v = v @ Y4 % p
        rows.append(v)
    return np.array(rows)
R = regwords()
print('rank of explicit regularised (u-derivative) basis at pair-by-pair ABJM point:', rank_mod(R))
pickle.dump({'R':R}, open('regL4.pkl','wb'))
# d-line valuation structure
coef = pickle.load(open('dline_coef.pkl','rb'))['coef']
N = 4
Ms = np.array(coef[:N+1])
sel, rest = valuation_basis(Ms)
print('d-line valuations', sorted(Counter([s[0] for s in sel]).items()), 'unresolved', len(rest))
# d-derivative regularised words: same recipe, with eps-Taylor coefficients of the words
Mw = {k: [coef[m][i] for m in range(N+1)] for i,k in enumerate(WIDX)}
def is_reg(a,k): return a == 3 or (a == 2 and k <= 2)
# per pair: the order-1 part of (w(a,k) - w(1,k+1)) or w(3,k); build products by working with full series rows
def dreg_row(a1,k2,a3,k4):
    # build combination series of words: product structure of per-pair differences
    terms1 = [((a1,k2),1)] if not (a1==2 and k2<=2) else [((2,k2),1),((1,k2+1),-1)]
    terms3 = [((a3,k4),1)] if not (a3==2 and k4<=2) else [((2,k4),1),((1,k4+1),-1)]
    n = (1 if is_reg(a1,k2) else 0) + (1 if is_reg(a3,k4) else 0)
    acc = np.zeros(256, dtype=np.int64)
    for (q1,c1) in terms1:
        for (q3,c3) in terms3:
            i = WIDX.index((q1[0],q1[1],q3[0],q3[1]))
            acc = (acc + c1*c3*coef[n][i]) % p
            # check lower orders cancel
    return acc, n
Rd = []; ns = []
for w in WIDX:
    r, n = dreg_row(*w); Rd.append(r); ns.append(n)
Rd = np.array(Rd)
print('rank of d-derivative regularised basis:', rank_mod(Rd), Counter(ns))
# compare d-derivative vs u-derivative rows: proportional?
prop = 0; inspan = 0
for i,w in enumerate(WIDX):
    if rank_mod(np.array([R[i], Rd[i]])) == 1: prop += 1
pickle.dump({'R':R,'Rd':Rd,'ns':ns}, open('regL4_both.pkl','wb'))
print('# rows where d-derivative row is proportional to u-derivative row:', prop, 'of 256')
