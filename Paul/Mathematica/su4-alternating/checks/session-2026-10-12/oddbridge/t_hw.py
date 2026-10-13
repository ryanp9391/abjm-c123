import numpy as np, time
from hw import *
Lp=3; N=6
d=hwdata()
S=d['L%d'%Lp]
secs=sectors(Lp)
print({k:len(v) for k,v in secs.items() if k in [(3,0,0,-3),(2,0,0,-2),(1,0,0,-1)]})
thA=np.array([1.,0,1,0,1,0]); z1=np.ones(4)
t0=time.time()
for i,st in enumerate(S):
    idx=secs[tuple(st['w'])]
    M=sector_matrix(Lp,idx,thA,z1)
    psi=st['psi'][idx]
    # ket eigenvalue test
    v=M@psi; lam=(psi.conj()@v)/(psi.conj()@psi)
    res=np.linalg.norm(v-lam*psi)
    # also HW check, norm outside sector
    print(i+1,st['dyn'],len(idx),'lam=%.6f%+.6fj'%(lam.real,lam.imag),'res %.2e'%res,'outside norm %.1e'%np.linalg.norm(np.delete(st['psi'],idx)))
print(time.time()-t0)
