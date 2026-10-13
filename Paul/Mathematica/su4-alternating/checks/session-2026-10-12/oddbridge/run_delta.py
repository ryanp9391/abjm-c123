import sys, time, pickle
sys.path.insert(0,'.')
from untw import *
np.set_printoptions(linewidth=200,precision=6)
def fmt(x):
    x=complex(x); return "%.6g%+.6gj"%(x.real,x.imag) if abs(x.imag)>1e-12 else "%.7g"%x.real
a=1/3; cc=2.0; phi=PHI1
IDX=[0,4,5,6]
for eps in [1e-2,1e-3,1e-4]:
    print("== eps=%g, family (a) pair 1 (k=2), phi1, c=2"%eps)
    for delta in [1.0,0.3,0.1,1e-2,1e-3,1e-4,1e-6,0.0]:
        th=THABJM+delta*THDIR; z=1+eps*phi
        tw=twisted(IDX,th,z)
        cells=[]
        for i in IDX:
            r=row(i,tw[i],2,a,cc,phi,eps,delta)
            psi=tw[i][0]
            c1=ccl(1,th); d2=dcl(2,th)
            t1_=tau(psi,th,z,th[0]); t2_=tau(psi,th,z,th[1])
            cells.append("st%d T=%s |ov|=%.6f 2nd=%.3f hat(tau1)=%s hat(tau2)=%s"%(i+1,fmt(r['T']),tw[i][1],tw[i][2],fmt(c1*t1_),fmt(t2_/d2)))
        print(" delta=%-7g c1=%.6g d2=%.6g | "%(delta,ccl(1,th),dcl(2,th))+" | ".join(cells))
