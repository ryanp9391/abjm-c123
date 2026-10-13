import numpy as np, sys
from ob import *
rng=np.random.default_rng(5)
def run(Lp, th, z, a=1/3):
    N=2*Lp
    sv0=Wv(a)*Lp
    Lam=np.diag(z); LamI=np.diag(1/np.array(z))
    # predecessor chain: P_k = first k sites rotated alternately
    def P(k):
        sv=list(sv0)
        for s in range(k): sv=rot(sv,s,LamI if s%2==0 else Lam)
        return prod_state(sv)
    out=[]
    for k in range(1,N+1):
        Pk=P(k); Pkm=P(k-1)
        if k%2==1:
            y=t1_apply(Pk, th, z, th[k-1])      # P_k t1(th_k) = (1/d_k) P_{k-1}
            r,dev,zc=ratio(y,Pkm); dk=1/r
        else:
            y=t1_apply(Pkm, th, z, th[k-1])     # P_{k-1} t1(th_k) = d_k P_k
            r,dev,zc=ratio(y,Pk); dk=r
        out.append((k,dk,d_closed(k,th,Lp),dev,zc))
    return out
for Lp in [1,2,3]:
    th=TH8 if Lp<4 else None
    z=Z4
    res=run(Lp, TH8[:2*Lp], Z4)
    print("L=",Lp)
    for k,dk,dc,dev,zc in res: print("  k=%d d=%.12g closed=%.12g relerr=%.2e prodDev=%.2e zeroComp=%.2e"%(k,dk,dc,abs(dk-dc)/abs(dc),dev,zc))
