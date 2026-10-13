import time, pickle
from taylor import *
from collections import Counter
F = ModP(); h = H; p = F.p
tF = Fr(-3,5); tA = tF + 2*h
th = [tA, tF, tA, tF]
ch = Chain(F, TYPES, th, h, G)
t0=time.time()
N = 6
TA = [t_taylor(ch, a, tA-h, N) for a in range(4)]
TF = t_taylor(ch, 1, tF, N)
print('taylor done', time.time()-t0, flush=True)
pickle.dump({'TA':TA,'TF':TF}, open('taylor_phys.pkl','wb'))
def run(lams, label):
    l1,l2,l3,l4 = lams
    ops = {'A1':[ser_scale_shift(TA[a], l1) for a in range(4)], 'A3':[ser_scale_shift(TA[a], l3) for a in range(4)],
           'T2': ser_scale_shift(TF, l2), 'T4': ser_scale_shift(TF, l4)}
    Ms = ser_words(ch, ops)
    sel, rest = valuation_basis(Ms)
    print(label, lams, 'valuations', sorted(Counter([s[0] for s in sel]).items()), 'unresolved', len(rest), time.time()-t0, flush=True)
    return sel
run((Fr(1),Fr(2),Fr(-3),Fr(5)), 'generic dirs')
run((Fr(0),Fr(0),Fr(1),Fr(1)), 'shift only pair 2 (A3,F4)')
run((Fr(0),Fr(1),Fr(0),Fr(-1)), 'shift only F sites')
