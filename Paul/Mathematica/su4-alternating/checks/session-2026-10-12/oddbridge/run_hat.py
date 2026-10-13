import sys
sys.path.insert(0,'.')
from untw import *
def fmt(x):
    x=complex(x); return "%.6g%+.6gj"%(x.real,x.imag) if abs(x.imag)>1e-12 else "%.7g"%x.real
a=1/3
PHI3=np.array([0.5,0.1,-0.4,-1.3])
IDX=[0,4,5,6]
for phn,phi in (("phi1",PHI1),("phi2",PHI2),("phi3",PHI3)):
    for eps in [1e-3,1e-4,1e-5]:
        th=THABJM; z=1+eps*phi   # delta=0
        tw=twisted(IDX,th,z)
        out=[]
        for i in IDX:
            psi=tw[i][0]
            c1=ccl(1,th); d2=dcl(2,th)
            ht1=c1*tau(psi,th,z,th[0]); ht2=tau(psi,th,z,th[1])/d2
            Wb=prod_state(Wv(a)*Lp).reshape(-1)
            ovW=Wb@psi
            K=ovW/(eps*a*(phi[3]-phi[0]))
            formA=(ht2-1)/ht1*K
            formB=(1-1/ht1)*K
            out.append("st%d hat1=%s hat2=%s K=%s  (t2-1)/t1*K=%s (1-1/t1)K=%s C1=%s"%(i+1,fmt(ht1),fmt(ht2),fmt(K),fmt(formA),fmt(formB),fmt(C1[i])))
        print(phn,"eps=%g delta=0"%eps); print("   "+"\n   ".join(out))
