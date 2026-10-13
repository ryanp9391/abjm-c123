import time, pickle
from taylor import *
F = ModP(); h = H; p = F.p
tF1, tF2 = Fr(-3,5), Fr(-11,8)
th = [tF1+2*h, tF1, tF2+2*h, tF2]
ch = Chain(F, TYPES, th, h, G)
t0=time.time()
N = 3
TA1 = [t_taylor(ch, a, th[0]-h, N) for a in range(4)]
TA3 = [t_taylor(ch, a, th[2]-h, N) for a in range(4)]
print('taylor done', time.time()-t0, flush=True)
pickle.dump({'TA1':TA1,'TA3':TA3}, open('taylor_pbp.pkl','wb'))
# checks: t3 coefficient 0 vanishes, t2_0 = t1_0 t1(thF)
Y2 = ch.t(1, th[1]); Y4 = ch.t(1, th[3])
print('t3_0 == 0', not np.any(TA1[3][0]), not np.any(TA3[3][0]))
print('t2_0 == t1_0 Y', not np.any((TA1[2][0]-F.mm(TA1[1][0],Y2))%p), not np.any((TA3[2][0]-F.mm(TA3[1][0],Y4))%p))
lam1, lam3 = Fr(1), Fr(1)
ops = {'A1':[ser_scale_shift(TA1[a], lam1) for a in range(4)], 'A3':[ser_scale_shift(TA3[a], lam3) for a in range(4)]}
T2 = np.zeros((N+1,256,256),dtype=np.int64); T2[0]=Y2
T4 = np.zeros((N+1,256,256),dtype=np.int64); T4[0]=Y4
ops['T2']=T2; ops['T4']=T4
Ms = ser_words(ch, ops)
print('series words', time.time()-t0, flush=True)
sel, rest = valuation_basis(Ms)
from collections import Counter
print('valuation counts', Counter([s[0] for s in sel]), 'unresolved', len(rest), time.time()-t0)
pickle.dump({'sel':[(s[0], s[1]) for s in sel]}, open('val_pbp.pkl','wb'))
