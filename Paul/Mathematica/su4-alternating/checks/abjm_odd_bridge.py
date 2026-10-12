import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'abjm_pair_base.py')).read().split('words={}')[0])
T1=[C4.t(1,th[s]) for s in range(4)]; inv=np.linalg.inv
X=inv(T1[0]); Y=T1[1]; Z=T1[2]; U=T1[3]
gl=prod([inv(Lam),Lam,inv(Lam),Lam])          # rho(gamma), gamma=Lam^-1 (A: Lam^-1, F: Lam)
Tw=X@Y@inv(Z)@U; c=np.trace(Tw@inv(gl))/256
print("t1(th1)^-1 t1(th2) t1(th3)^-1 t1(th4) = c rho(gamma) as OPERATOR:",np.allclose(Tw,c*gl), " c=%.6g"%c)
par=lambda x,y: abs(abs(np.vdot(x,y))-np.linalg.norm(x)*np.linalg.norm(y))<1e-9*np.linalg.norm(x)*np.linalg.norm(y)
PL=X@Y; PR=Z@inv(U)
print("<W|PL PR = rho_12(gam) rho_34(gam^-1) <W| ?",par(W@PL@PR,prod([vA@inv(Lam),wF@Lam,vA@Lam,wF@inv(Lam)])))
print("<W|PL^2 product?",par(W@PL@PL,prod([vA@inv(Lam)@inv(Lam),wF@Lam@Lam,vA,wF])))
# odd bridge l=1: site-1 kappa family, kappa = z4 gives vacuum A field on site 1
t=lambda i,j: th[i]-th[j]; c1=1/(h*t(0,1)*(t(0,2)+h)*t(0,3))
print("<W| t1(th1)^-1 = c1 <W|rho_1(Lam^-1), c1=1/(h th12 (th13+h) th14):",np.allclose(W@X,c1*prod([vA@inv(Lam),wF,vA,wF])))
kap=z[3]; K=np.eye(256)-kap*X/c1
tg=prod([vA@(np.eye(4)-kap*inv(Lam)),wF,vA,wF])
print("<W|(1 - z4 t1(th1)^-1/c1) =",np.round((vA@(np.eye(4)-kap*inv(Lam))),4),"(x) W-rest :",np.allclose(W@K,tg))
