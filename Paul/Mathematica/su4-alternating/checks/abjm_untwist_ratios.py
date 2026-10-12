import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'abjm_untwist.py')).read().split('for eps in')[0])
def run(eps,phi):
    z=1+eps*phi; C=chain(z); V=eigs(C); T1=C.t(1,th[0]); T2=C.t(1,th[1])
    b=z[0]/(z[0]-z[3]); ref=Wst([b-1,b-1]); loc=Wst([b,b-1])
    c=np.vdot(loc,ref@inv(T1)@T2)/np.vdot(loc,loc)
    out=[]
    for k in range(256):
        v=V[:,k]; v=v/np.linalg.norm(v)
        r=max(np.linalg.norm(Jm@v) for Jm in Jr)
        if r>30*eps: continue
        T=tgt@v
        if abs(T)<1e-6: continue
        F=(1/c)*(ref@inv(T1)@T2@v)
        # identify the HW limit by its untwisted t1 eigenvalue at a test point
        out.append((round(float(np.real((v.conj()@C0.t(1,0.123)@v))),4), F/T))
    return sorted(out)
for phi_ in (np.array([0.9,0.4,-0.3,-1.1]),np.array([1.3,-0.2,0.5,-0.7])):
    print("phi=",phi_)
    for eps in (1e-3,1e-4):
        print("  eps=%.0e:"%eps,["%.5f%+.5fi"%(r.real,r.imag) for _,r in run(eps,phi_)])
