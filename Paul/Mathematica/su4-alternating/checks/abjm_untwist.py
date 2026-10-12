import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'rotbase.py')).read().split('kap=0.37')[0])
inv=np.linalg.inv
pairW=lambda a:(e1+a*e4, e4-a*e1)
Wst=lambda alist: prod([x for a in alist for x in pairW(a)])
vac2=prod([e1,e4,e1,e4]); tgt=prod([e1,e4,*pairW(-1.0)])
# raising operators on kets: J_ij = sum_s rho_s(E_ij)
reps=[rA,rF,rA,rF]
def J(i,j): return sum(functools.reduce(np.kron,[reps[s](i,j) if b==s else I for b in range(4)]) for s in range(4))
Jr=[J(i,j) for i in range(4) for j in range(i+1,4)]
def chain(zv):
    global G; G=np.diag(zv); C=Chain([rA,rF,rA,rF],th)
    return C
def eigs(C):
    M=C.t(1,0.3711)+0.37*C.t(2,-0.917)+0.11*C.t(3,0.523)
    w,V=np.linalg.eig(M); return V
phi=np.array([0.9,0.4,-0.3,-1.1])
# untwisted HW check of b-independence
C0=chain(np.ones(4)); V0=eigs(C0)
hw=[k for k in range(256) if max(np.linalg.norm(Jm@V0[:,k]) for Jm in Jr)<1e-8*np.linalg.norm(V0[:,k])]
print("untwisted: eigvecs that are HW:",len(hw),"(degenerate multiplets make eig basis mixed; using projector instead)")
# HW subspace = common kernel of raising ops, restricted per eigenvalue: build kernel directly
K=np.vstack(Jr); _,s,Vh=np.linalg.svd(K); HW=Vh[np.sum(s>1e-9):].T; print("dim HW subspace:",HW.shape[1])
for b in (0.7,-2.3,5.0):
    lhs=tgt@HW; rhs=Wst([b,b-1])@HW; print("  b=%.1f  max|<vac x W_-1|hw> - <W_b x W_{b-1}|hw>| = %.1e"%(b,np.abs(lhs-rhs).max()))
for eps in (1e-1,1e-2,1e-3,1e-4):
    z=1+eps*phi; C=chain(z); V=eigs(C)
    T1=C.t(1,th[0]); T2=C.t(1,th[1])
    b=z[0]/(z[0]-z[3])
    ref=Wst([b-1,b-1]); loc=Wst([b,b-1])
    # exactness of the localisation identity at finite twist
    lhsop=ref@inv(T1)@T2; c=np.vdot(loc,lhsop)/np.vdot(loc,loc)
    ok=np.linalg.norm(lhsop-c*loc)/np.linalg.norm(lhsop)
    # states close to HW
    hwlike=[k for k in range(256) if max(np.linalg.norm(Jm@V[:,k]) for Jm in Jr)<30*eps*np.linalg.norm(V[:,k])]
    ratios=[]
    for k in hwlike:
        v=V[:,k]; T=tgt@v
        if abs(T)<1e-6*np.linalg.norm(v): continue
        tau1=(vac2@T1@v)/(vac2@v) if abs(vac2@v)>1e-12 else None
        F=(1/c)*(ref@inv(T1)@T2@v)
        ratios.append(F/T)
    ratios=np.array(ratios)
    print("eps=%.0e  b=%.3g  localisation resid %.0e  HW-like states with T!=0: %d   |F/T - 1| max %.2e median %.2e"%(eps,b,ok,len(ratios),np.abs(ratios-1).max(),np.median(np.abs(ratios-1))))
