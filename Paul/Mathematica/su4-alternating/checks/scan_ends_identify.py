import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'scan_classify.py')).read().split('want=[')[0])
k1,k2=0.37,-0.81; v=e1@(I-k1*G.T)@(I-k2*G.T); refs=[v,e4,v,e4]; Om=prod(refs)
gam=np.linalg.inv(G).T; gi=G.T          # gamma=G^{-T}; on A acts by gamma, on F by gamma^{-T}=G
A=lambda g: v@g; F=lambda g: e4@np.linalg.inv(g).T
def par(x,y): return np.linalg.matrix_rank(np.array([x,y]),1e-9*np.linalg.norm(x)*0+1e-9)==1 if False else abs(abs(np.vdot(x,y))-np.linalg.norm(x)*np.linalg.norm(y))<1e-9*np.linalg.norm(x)*np.linalg.norm(y)
o=ops; X=o['t1(th1+0 h)^-1']; Y=o['t1(th2+0 h)']; Zi=o['t1(th3+0 h)^-1']; Z=o['t1(th3+0 h)']; Yi=o['t1(th2+0 h)^-1']; Wi=o['t1(th4+0 h)^-1']; Wp=o['t1(th4+0 h)']
tests={
 "X^-1         = rho_1(gam)":          (Om@X,        prod([A(gam),e4,v,e4])),
 "X^-1 Y       = rho_12(gam)":         (Om@X@Y,      prod([A(gam),F(gam),v,e4])),
 "X^-1 Y Z^-1  = rho_123(gam)":        (Om@X@Y@Zi,   prod([A(gam),F(gam),A(gam),e4])),
 "X^-1 Y Z^-1 W = rho(gam) global":    (Om@X@Y@Zi@Wp,prod([A(gam),F(gam),A(gam),F(gam)])),
 "W^-1         = rho_4(gam^-1)":       (Om@Wi,       prod([v,e4,v,F(gi)])),
 "Z W^-1       = rho_34(gam^-1)":      (Om@Z@Wi,     prod([v,e4,A(gi),F(gi)])),
 "Y^-1 Z W^-1  = rho_234(gam^-1)":     (Om@Yi@Z@Wi,  prod([v,F(gi),A(gi),F(gi)])),
}
for k,(x,y) in tests.items(): print("%-36s %s"%(k,par(x,y)))
