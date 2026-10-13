import sys, pickle
sys.path.insert(0,'.')
from untw import *
def f(x):
    x=complex(x)
    if abs(x)<1e-13: return "0"
    if abs(x.imag)<1e-13*max(1,abs(x.real)): return "%.6g"%x.real
    if abs(x.real)<1e-13*max(1,abs(x.imag)): return "%.6gi"%x.imag
    return "%.5g%+.5gi"%(x.real,x.imag)
a=1/3
dyn=lambda i: "[%d,%d,%d]"%tuple(S[i]['dyn'])
ids=[0,4,5,6]
out=[]
out.append("#### delta-scan: F-site family, pair 1, phi1, c = 2, eps = 1e-3 (theta = theta_ABJM + delta theta0)")
out.append("| delta | state (E) | c1 | d2 | |<Psi0|Psi_eps>| | hat tau_1 | hat tau_2 | T_eps | limit | |T - limit| |")
out.append("|---|---|---|---|---|---|---|---|---|---|")
eps=1e-3
for dl in [0.1,0.01,1e-3,1e-4,1e-6,0.0]:
    th=THABJM+dl*THDIR; z=1+eps*PHI1
    tw=twisted(ids,th,z)
    for i in ids:
        psi=tw[i][0]; taus=[tau(psi,th,z,th[j]) for j in range(2)]
        r=row(i,tw[i],2,a,2.0,PHI1,eps,dl)
        c1=ccl(1,th); d2=dcl(2,th)
        out.append("| %g | %s (%g) | %.6g | %.6g | %.6f | %s | %s | %s | %s | %.2e |"%(dl,dyn(i),S[i]['E'],c1,d2,tw[i][1],f(c1*taus[0]),f(taus[1]/d2),f(r['T']),f(r['lim']),abs(r['T']-r['lim'])))
SC="\n".join(out)
out=[]
out.append("#### delta = 0 (exactly homogeneous, twisted): hat tau_1, hat tau_2, K = <W|Psi_eps>/(eps a (phi4 - phi1)) and the c-free formula C_1 = (1 - 1/hat tau_1) K")
out.append("| phi | eps | state (E) | hat tau_1 | hat tau_2 | K | (1 - 1/hat tau_1) K | (hat tau_2 - 1) K/hat tau_1 | C_1 (data) |")
out.append("|---|---|---|---|---|---|---|---|---|")
PHI3=np.array([0.5,0.1,-0.4,-1.3])
Wb=prod_state(Wv(a)*Lp).reshape(-1)
th=THABJM
for pn,phi in (("phi1",PHI1),("phi2",PHI2),("phi3",PHI3)):
    for eps in (1e-3,1e-4,1e-5):
        z=1+eps*phi
        tw=twisted([4,5,6],th,z)
        for i in (4,5,6):
            psi=tw[i][0]; t1=ccl(1,th)*tau(psi,th,z,th[0]); t2=tau(psi,th,z,th[1])/dcl(2,th)
            K=(Wb@psi)/(eps*a*(phi[3]-phi[0]))
            out.append("| %s | %g | %s (%g) | %s | %s | %s | %s | %s | %s |"%(pn,eps,dyn(i),S[i]['E'],f(t1),f(t2),f(K),f((1-1/t1)*K),f((t2-1)/t1*K),f(C1[i])))
HT="\n".join(out)
# vacuum
out=[]
out.append("#### twisted vacuum, L = 2, theta = (0.3, 0.3-d, 0.3, 0.3-d), z = (1.3, 0.8, 1.1, 0.7): tau(theta_k), d_k, hat tau_k")
out.append("| d = theta_A - theta_F | k | tau(theta_k) on vac | d_k | hat tau_k | z_1 (k odd) or z_4 (k even) |")
out.append("|---|---|---|---|---|---|")
import ob
z=np.array([1.3,0.8,1.1,0.7]); vac=ob.prod_state([ob.E(1),ob.E(4)]*2)
for d in (1.0,0.1,1e-3,1e-5):
    thv=np.array([0.3,0.3-d,0.3,0.3-d])
    for k in range(1,5):
        y=ob.t1_apply(vac,thv,z,thv[k-1]); lam=(y.reshape(-1)@vac.reshape(-1))
        dk=ob.d_closed(k,thv,2); hat=dk*lam if k%2==1 else lam/dk
        out.append("| %g | %d | %.5g | %.5g | %.10g | %g |"%(d,k,lam,dk,hat,z[0] if k%2==1 else z[3]))
VT="\n".join(out)
pickle.dump(dict(SC=SC,HT=HT,VT=VT),open('tables_part3.pkl','wb'))
print(SC); print(HT); print(VT)
