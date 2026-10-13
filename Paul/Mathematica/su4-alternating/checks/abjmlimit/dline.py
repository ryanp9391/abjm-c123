import time, pickle, sys
from taylor import *
F = ModP(); h = H; p = F.p
tF1, tF2 = Fr(-3,5), Fr(-11,8)
def th_of(e): return [tF1+2*h+e, tF1, tF2+2*h+e, tF2]
npts = 27
nodes = [Fr(s) for s in range(1, npts+1)]
t0 = time.time()
vals = []
for x in nodes:
    ch = Chain(F, TYPES, th_of(x), h, G)
    vals.append(words(ch))
    del ch
print('evaluated', time.time()-t0, flush=True)
# interpolate: coefficients c_m = sum_i Vinv[m][i] vals[i]
from taylor import frac_inv
V = [[x**k for k in range(npts)] for x in nodes]
Vi = frac_inv(V)
coef = []
for m in range(npts):
    acc = np.zeros((256,256), dtype=np.int64)
    for i in range(npts):
        w = F.c(Vi[m][i])
        if w: acc = (acc + w*vals[i]) % p
    coef.append(acc)
print('degrees: nonzero top coeffs', [m for m in range(npts) if np.any(coef[m])][-3:], time.time()-t0)
# check consistency with direct evaluation at eps=0
ch0 = Chain(F, TYPES, th_of(Fr(0)), h, G)
M0 = words(ch0)
print('coef0 == words at eps=0:', not np.any((coef[0]-M0)%p))
pickle.dump({'coef': coef[:8]}, open('dline_coef.pkl','wb'))
