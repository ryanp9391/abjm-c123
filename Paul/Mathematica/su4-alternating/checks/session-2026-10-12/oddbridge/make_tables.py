import sys, pickle
sys.path.insert(0,'.')
from untw import *
def f(x):
    x=complex(x)
    if abs(x)<1e-13: return "0"
    if abs(x.imag)<1e-13*max(1,abs(x.real)): return "%.5g"%x.real
    if abs(x.real)<1e-13*max(1,abs(x.imag)): return "%.5gi"%x.imag
    return "%.4g%+.4gi"%(x.real,x.imag)
out=[]
a=1/3
ids=[0,1,2,4,5,6]       # WL ids 1,2,3,5,6,7
order=[0,4,5,6,1,2]     # obStates = {1,5,6,7,2,3}
dyn=lambda i: "[%d,%d,%d]"%tuple(S[i]['dyn'])
# --- state table
out.append("### STATE TABLE")
out.append("| # | Dynkin | weight | sector dim | E_MZ | V0 | X11_p (p=1,2,3) | X44_p | C_1 (data) |")
out.append("|---|---|---|---|---|---|---|---|---|")
for i,s in enumerate(S):
    out.append("| %d | %s | %s | %d | %.4g | %s | %s | %s | %s |"%(i+1,dyn(i),tuple(int(x) for x in s['w']),len(SECS[tuple(s['w'])]),s['E'],f(V0[i]),", ".join(f(x) for x in X11[i]),", ".join(f(x) for x in X44[i]),f(C1[i])))
# --- twisted & rows
EPS=[1e-2,1e-3,1e-4]
TW={}
for phn,phi in (("phi1",PHI1),("phi2",PHI2)):
    for eps in EPS:
        th=THABJM+eps*THDIR; z=1+eps*phi
        tw=twisted(ids,th,z)
        for i in ids:
            tw[i]=list(tw[i])+[ [tau(tw[i][0],th,z,th[j]) for j in range(6)] ]
        TW[(phn,eps)]=tw
pickle.dump(TW,open('tw2.pkl','wb'))
def mkrow(i,k,cc,phi,phn,eps):
    th=THABJM+eps*THDIR; z=1+eps*phi
    tw=TW[(phn,eps)][i]; r=row(i,tw[:4],k,a,cc,phi,eps,eps)
    return r
out.append("### A1/A3 table (phi1, c=2, a=1/3)")
out.append("| k (pair) | family | Dynkin | E_MZ | eps | ov | ov2 | <W|Psi_eps> | T_eps | limit | err | idres |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for k in (2,1,3,4,5,6):
    for eps in EPS:
        for i in order:
            r=mkrow(i,k,2.0,PHI1,"phi1",eps); tw=TW[("phi1",eps)][i]
            out.append("| %d (%d) | %s | %s | %.4g | %g | %.6f | %.4f | %s | %s | %s | %.2e | %.1e |"%(k,(k+1)//2,"F (a)" if k%2==0 else "A (b)",dyn(i),S[i]['E'],eps,tw[1],tw[2],f(r['ovW']),f(r['T']),f(r['lim']),abs(r['T']-r['lim']),abs(r['lhs_prod']-r['rhs_tau'])))
pickle.dump(out,open('tables_part1.pkl','wb'))
open('tables_part1.txt','w',encoding='utf8').write("\n".join(out))
print("ok")
