import time, sys
from words import *
F = ModP()
h = H
def rk(th):
    ch = Chain(F, TYPES, th, h, G)
    return rank_mod(words(ch))
tF1, tF2 = Fr(-3,5), Fr(-11,8)
out = []
ds = [Fr(1,3), Fr(0), h, 2*h, 3*h, 4*h, 2*h+Fr(1,1000), 2*h-Fr(1,1000), 2*h+Fr(1,10**6)]
for d in ds:
    r1 = rk([tF1+d, tF1, tF2+d, tF2])
    r2 = rk([tF1+d, tF1, tF1+d, tF1])
    print('d=',d, ' pair-by-pair rank', r1, ' pair-periodic rank', r2, flush=True)
