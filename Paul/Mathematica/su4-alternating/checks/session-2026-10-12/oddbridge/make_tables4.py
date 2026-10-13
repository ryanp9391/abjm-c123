import sys, pickle, re, math
sys.path.insert(0,'.')
from untw import *
def f(x):
    x=complex(x)
    if abs(x)<1e-13: return "0"
    if abs(x.imag)<1e-13*max(1,abs(x.real)): return "%.5g"%x.real
    if abs(x.real)<1e-13*max(1,abs(x.imag)): return "%.5gi"%x.imag
    return "%.4g%+.4gi"%(x.real,x.imag)
a=1/3
TW=pickle.load(open('tw2.pkl','rb'))
dyn=lambda i: "[%d,%d,%d]"%tuple(S[i]['dyn'])
def mkrow(i,k,cc,phi,phn,eps):
    tw=TW[(phn,eps)][i]; return row(i,tw[:4],k,a,cc,phi,eps,eps)
out=["| k (pair) | family | state (E_MZ) | eps | T_eps | limit | |T - limit| | |LHS_prod - RHS_tau| |","|---|---|---|---|---|---|---|---|"]
for k in (3,4,5,6):
    for i in (0,4,5,6):
        for eps in (1e-3,1e-4):
            r=mkrow(i,k,2.0,PHI1,"phi1",eps)
            out.append("| %d (%d) | %s | %s (%g) | %g | %s | %s | %.2e | %.1e |"%(k,(k+1)//2,"F (a)" if k%2==0 else "A (b)",dyn(i),S[i]['E'],eps,f(r['T']),f(r['lim']),abs(r['T']-r['lim']),abs(r['lhs_prod']-r['rhs_tau'])))
print("\n".join(out))
print()
# rates for k=3..6
out=["| k | state | rate 1e-2 -> 1e-3 | rate 1e-3 -> 1e-4 |","|---|---|---|---|"]
for k in (3,4,5,6):
    for i in (0,4,6):
        e=[]
        for eps in (1e-2,1e-3,1e-4):
            r=mkrow(i,k,2.0,PHI1,"phi1",eps); e.append(abs(r['T']-r['lim']))
        out.append("| %d | %s (%g) | %.3f | %.3f |"%(k,dyn(i),S[i]['E'],math.log10(e[0]/e[1]),math.log10(e[1]/e[2])))
print("\n".join(out))
