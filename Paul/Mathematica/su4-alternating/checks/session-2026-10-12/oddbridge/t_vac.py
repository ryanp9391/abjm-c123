import sys; sys.path.insert(0,'.')
from ob import *
# eigenvalue of tilde t_k = d_k^{+-1} t1(th_k) on the (twisted) vacuum bra vac^L; theta_A - theta_F = d -> 0
for Lp in (2,3):
    z=np.array([1.3,0.8,1.1,0.7])
    vac=prod_state([E(1),E(4)]*Lp)
    for d in (1.0,0.5,0.1,1e-2,1e-3,1e-4,0.0):
        thA=np.array([0.3+0.11*i for i in range(Lp)]); 
        th=[]
        for i in range(Lp): th+= [thA[i], thA[i]-d + 0.07*i*(d!=0.0)]
        th=np.array(th)
        out=[]
        for k in (1,2,3):
            y=t1_apply(vac,th,z,th[k-1]); lam=ratio(y,vac)[0]
            dk=d_closed(k,th,Lp) if abs(d)>0 or True else 0
            out.append("t%d=%.6g d%d=%.4g"%(k,lam,k,dk) if abs(dk)<1e12 else "t%d=%.6g d%d=inf"%(k,lam,k))
        try:
            hat=[ratio(t1_apply(vac,th,z,th[0]),vac)[0]*d_closed(1,th,Lp)]
        except ZeroDivisionError: hat=['inf']
        print(Lp,"d=%g"%d,"|"," ".join(out),"| hat1=",hat)
