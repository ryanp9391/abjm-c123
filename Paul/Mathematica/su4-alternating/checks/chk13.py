import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'chk7.py')).read().split('r,c=prop(vac@Y')[0])
c1=1/(h*t(0,1)*(t(0,2)+h)*t(0,3)); l3=h*t(2,1)*(t(2,0)+h)*t(2,3)
o1=X/c1; o3=Z/l3
P13=vac@so(0,gam)@so(2,G.T)              # sites 1 and 3 both rotated
word=vac@((t(0,2)-h)*o1@o3+h*np.eye(256))/t(0,2)
print("P13 = [(th13-h) o1 o3 + h]/th13 :  resid %.1e"%(np.linalg.norm(P13-word)/np.linalg.norm(P13)))
v=vac@o1@o3; print("<0|o1 o3 is a pure product rotation? resid %.2f"%(np.linalg.norm(v-P13*np.vdot(P13,v)/np.vdot(P13,P13))/np.linalg.norm(v)))
# eigenstate check of overlap formula
w,V=np.linalg.eig(C4.t(1,0.123)+0.7*C4.t(2,-0.31)); k1,k3=0.41,-0.83
cov=functools.reduce(np.kron,[e1@(I-k1*gam),e4,e1@(I-k3*G.T),e4])
o1h=np.array([ (vac@o1@V[:,i])/(vac@V[:,i]) for i in range(256)]); o3h=np.array([(vac@o3@V[:,i])/(vac@V[:,i]) for i in range(256)])
pred=((1-k1*o1h)*(1-k3*o3h)-(h/t(0,2))*k1*k3*(o1h*o3h-1))*(vac@V)
print("overlap formula on 256 eigenstates: max rel err %.1e"%np.max(np.abs(cov@V-pred)/np.abs(cov@V)))
err=np.abs(cov@V-pred)/np.abs(cov@V); res=np.array([np.linalg.norm((C4.t(1,0.123)+0.7*C4.t(2,-0.31))@V[:,i]-w[i]*V[:,i]) for i in range(256)])
print("median rel err %.1e; worst state: |<0|Psi>| %.1e, eigvec residual %.1e"%(np.median(err),abs(vac@V[:,err.argmax()]),res[err.argmax()]))
# direct operator-level version: cov = <0|[(1-k1 o1)(1-k3 o3) - (h/th13)k1k3(o1 o3 - 1)]
op=(np.eye(256)-k1*o1)@(np.eye(256)-k3*o3)-(h/t(0,2))*k1*k3*(o1@o3-np.eye(256))
print("covector identity: resid %.1e"%(np.linalg.norm(vac@op-cov)/np.linalg.norm(cov)))
