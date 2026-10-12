import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'rotbase.py')).read().split('kap=0.37')[0])
z=np.array([1.5,1.25,2/3,4/7]); Gc=G.copy(); Lam=np.diag(z)
V=np.array([[zj**(3-i) for zj in z] for i in range(4)])   # columns = right eigvecs (z^3,z^2,z,1)
print("Gc V = V Lam:",np.allclose(Gc@V,V@Lam))
u=0.3719
tc=C4.t(1,u); G=Lam; tL=C4.t(1,u); G=Gc
def RR(M): return prod_m([M,np.linalg.inv(M).T,M,np.linalg.inv(M).T])
prod_m=lambda Ms: functools.reduce(np.kron,Ms)
for nm,M in (("V",V),("V^-1",np.linalg.inv(V)),("V^T",V.T),("V^-T",np.linalg.inv(V).T)):
    R=RR(M); print("t^Gc R = R t^Lam with M=%-5s: %s"%(nm,np.allclose(tc@R,R@tL)))
