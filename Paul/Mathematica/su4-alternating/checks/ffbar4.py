import os
import numpy as np, functools, itertools
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ffbar2.py')).read().split('# L=2: diagonal')[0])
def expm(A):
    w,V=np.linalg.eig(A); return (V@np.diag(np.exp(w))@np.linalg.inv(V)).real
th=[0.2857,-0.6,0.731,-1.37]; reps=[rf,rb,rf,rb]
t1,so=chain(reps,th); I=np.eye(N)
rng=np.random.default_rng(3); g=expm(0.3*rng.standard_normal((N,N)))
gb=np.linalg.inv(g).T
rg=functools.reduce(np.kron,[g,gb,g,gb])
e1=np.eye(N)[0]; e4=np.eye(N)[3]
vac=functools.reduce(np.kron,[e1,e4,e1,e4]); Om=vac@rg
Ti={a:np.linalg.inv(t1(th[a])) for a in range(4)}
def prop(x,y):
    c=np.vdot(y,x)/np.vdot(y,y); return np.linalg.norm(x-c*y)/np.linalg.norm(x)
Gi=np.linalg.inv(G)
loc={0:so(0,Gi.T),1:so(1,Gi),2:so(2,Gi.T),3:so(3,Gi)}
for a in range(4):
    print("site",a+1,"<Om|t(th)^-1 vs <Om|rho_a(twist^-1): resid %.2e"%prop(Om@Ti[a],Om@loc[a]))
# words of length 2 vs block rotations
for a,b in itertools.permutations(range(4),2):
    r=prop(Om@Ti[a]@Ti[b],Om@loc[a]@loc[b])
    if r<1e-8: print("block",(a+1,b+1),"OK")
# rank of span of <Om| words (length<=2, inverses) projected: how many site-local covectors?
