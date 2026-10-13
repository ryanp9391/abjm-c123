import time
from taylor import *
F = ModP(); h = H; p = F.p
tF1, tF2 = Fr(-3,5), Fr(-11,8)
th = [tF1+2*h, tF1, tF2+2*h, tF2]
ch = Chain(F, TYPES, th, h, G)
t0=time.time()
# generic shifts at fixed theta
for sh in [(Fr(1,7), Fr(-1,5)), (Fr(-3*1,2), Fr(-3,2)), (Fr(1,1), Fr(1,1))]:
    M = words(ch, uA1=th[0]-h+sh[0], uA3=th[2]-h+sh[1])
    print('shift', sh, 'rank', rank_mod(M), time.time()-t0, flush=True)
