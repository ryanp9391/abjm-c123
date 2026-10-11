import os
import numpy as np, functools, itertools
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ffbar4.py')).read().split('Ti={a')[0])
Gi=np.linalg.inv(G)
def prop(x,y):
    c=np.vdot(y,x)/np.vdot(y,y); return np.linalg.norm(x-c*y)/np.linalg.norm(x), c
r,c=prop(Om@np.linalg.inv(t1(th[0]))@t1(th[1]), Om@so(0,Gi.T)@so(1,G))
t=lambda a,b: th[a]-th[b]
print("end block {1,2}: resid %.1e  c=%.8f  predicted %.8f"%(r,c.real, t(1,2)*(t(1,3)-h)/((t(0,2)+h)*t(0,3))))
r,c=prop(Om@np.linalg.inv(t1(th[3]))@t1(th[2]), Om@so(2,G.T)@so(3,Gi))
print("end block {3,4}: resid %.1e  c=%.8f  predicted %.8f"%(r,c.real, t(2,1)*(t(2,0)+h)/((t(3,1)-h)*t(3,0))))
