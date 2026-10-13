import time, pickle, sys
from taylor import *
from collections import Counter
F = ModP(); h = H; p = F.p
tF = Fr(-3,5)
dval = Fr(sys.argv[1])
N = int(sys.argv[2])
tA = tF + dval
th = [tA, tF, tA, tF]
ch = Chain(F, TYPES, th, h, G)
t0=time.time()
TA = [t_taylor(ch, a, tA-h, N) for a in range(4)]
TF = t_taylor(ch, 1, tF, N)
res = {}
for lams in [(Fr(1),Fr(2),Fr(-3),Fr(5)), (Fr(2),Fr(-1),Fr(7),Fr(3))]:
    l1,l2,l3,l4 = lams
    ops = {'A1':[ser_scale_shift(TA[a], l1) for a in range(4)], 'A3':[ser_scale_shift(TA[a], l3) for a in range(4)],
           'T2': ser_scale_shift(TF, l2), 'T4': ser_scale_shift(TF, l4)}
    Ms = ser_words(ch, ops)
    sel, rest = valuation_basis(Ms)
    res[lams] = (sorted(Counter([s[0] for s in sel]).items()), len(rest))
    print('d =', dval, 'directions', [str(x) for x in lams], 'valuation counts', res[lams][0], 'unresolved', len(rest), 'sum n_j =', sum(k*v for k,v in res[lams][0]), round(time.time()-t0,1), flush=True)
