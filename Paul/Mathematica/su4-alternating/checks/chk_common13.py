import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'chk7.py')).read().split('r,c=prop(vac@Y')[0])
T1=C4.t(1,th[0]); T3=C4.t(1,th[2])
w0=vac; w1=vac@T1; w3=vac@T3; w13=vac@T1@T3
R1=vac@so(0,G.T); R3=vac@so(2,G.T); R13=vac@so(0,G.T)@so(2,G.T)
B=np.array([R1,R3,R13]).T
def coords(w):
    rest=w-w0*(np.vdot(w0,w)/np.vdot(w0,w0)) if False else w
    A=np.array([w0,R1,R3,R13]).T; c,*_=np.linalg.lstsq(A,w,rcond=None)
    return c, np.linalg.norm(A@c-w)/np.linalg.norm(w)
for n,w in (("t1(th1)",w1),("t1(th3)",w3),("t1(th1)t1(th3)",w13)):
    c,r=coords(w); print("<0|%-15s = %s . (vac, R1, R3, R13)   resid %.1e"%(n,np.round(c,6),r))
# solve 1 + a t1(th1) + b t1(th3) + c t1t3 = lam * <0| rho1(g) rho3(g), g = 1 - kappa G^T
for kappa in (0.3,-1.7,2.4):
    target=vac-kappa*(R1+R3)+kappa**2*R13
    A=np.array([w1,w3,w13,-target]).T; sol,*_=np.linalg.lstsq(A,-w0,rcond=None)
    a,b,c,lam=sol; res=np.linalg.norm(w0+a*w1+b*w3+c*w13-lam*target)/np.linalg.norm(target)
    g=I-kappa*G.T; prod=functools.reduce(np.kron,[e1@g,e4,e1@g,e4])
    print("kappa=%5.2f  a=%.6g b=%.6g c=%.6g  resid %.1e  | target == <0|rho1(g)rho3(g): %.1e"%(kappa,a,b,c,res,np.linalg.norm(target-prod)))
