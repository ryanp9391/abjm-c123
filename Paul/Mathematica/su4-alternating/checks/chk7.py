import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'g2test.py')).read().split('rng=np.random')[0])
th=[2/7,-3/5,8/11,-11/8]; C4=Chain([rA,rF,rA,rF],th)
e1=I[0]; e4=I[3]; vac=functools.reduce(np.kron,[e1,e4,e1,e4])
Gi=np.linalg.inv(G); gam=Gi.T
def prop(x,y):
    c=np.vdot(y,x)/np.vdot(y,y); return np.linalg.norm(x-c*y)/np.linalg.norm(x),c.real
t=lambda a,b: th[a]-th[b]
X=np.linalg.inv(C4.t(1,th[0])); Y=C4.t(1,th[1]); Z=C4.t(1,th[2]); W=np.linalg.inv(C4.t(1,th[3]))
so=C4.so
r,c=prop(vac@Y,vac@so(1,G)); print("site2 alone: resid %.1e c=%.6g  l2=%.6g"%(r,c,h*t(0,1)*t(1,2)*(t(1,3)-h)))
r,c=prop(vac@Z,vac@so(2,G.T)); print("site3 alone: resid %.1e c=%.6g  l3=%.6g"%(r,c,h*t(2,1)*(t(2,0)+h)*t(2,3)))
c1=1/(h*t(0,1)*(t(0,2)+h)*t(0,3)); l2=h*t(0,1)*t(1,2)*(t(1,3)-h)
k1,k2=0.37,-1.21
lhs=vac@(np.eye(256)-k1*X/c1)@(np.eye(256)-k2*Y/l2)
v1=e1@(I-k1*gam); v2=e4@(I-k2*G)
rhs=functools.reduce(np.kron,[v1,v2,e1,e4]); print("block (k1,k2) family: resid %.1e"%(np.linalg.norm(lhs-rhs)/np.linalg.norm(rhs)))
