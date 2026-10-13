import sys, pickle
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
EPS=[1e-2,1e-3,1e-4]
def mkrow(i,k,cc,phi,phn,eps):
    tw=TW[(phn,eps)][i]; return row(i,tw[:4],k,a,cc,phi,eps,eps)
out=[]
# P1: pair-1 families compact
out.append("#### P1: pair 1 (k = 2 F-site family (a); k = 1 A-site family (b)), phi1, c = 2, a = 1/3")
out.append("| family | state (Dynkin, E_MZ) | eps | |<Psi0|Psi_eps>| | <W_a^3|Psi_eps> | T_eps | limit | |T - limit| | |LHS_prod - RHS_tau| |")
out.append("|---|---|---|---|---|---|---|---|---|")
for k in (2,1):
    for i in (0,4,5,6):
        for eps in EPS:
            r=mkrow(i,k,2.0,PHI1,"phi1",eps); tw=TW[("phi1",eps)][i]
            out.append("| %s | %s, %g | %g | %.6f | %s | %s | %s | %.2e | %.1e |"%("F (a), k=2" if k==2 else "A (b), k=1",dyn(i),S[i]['E'],eps,tw[1],f(r['ovW']),f(r['T']),f(r['lim']),abs(r['T']-r['lim']),abs(r['lhs_prod']-r['rhs_tau'])))
out.append("Controls (states 2 = [4,0,0] and 3 = [1,2,1], weight sectors (3,-1,-1,-1) and (2,1,-1,-2)): T_eps = 0 and limit = 0 identically (the product bras only have weights (x,0,0,-x)); the 18 remaining saved states (sectors up to 256) were checked in numpy: T_eps -> 0 (exactly 0 or O(eps) / O(eps^2)), identity residual <= 4e-15.")
P1="\n".join(out)
# rates
out=[]
out.append("#### rates log10(err(eps)/err(eps/10)) from the P1 rows (1 = O(eps))")
out.append("| family | state | rate 1e-2 -> 1e-3 | rate 1e-3 -> 1e-4 |")
out.append("|---|---|---|---|")
import math
for k in (2,1):
    for i in (0,4,5,6):
        e=[abs(mkrow(i,k,2.0,PHI1,"phi1",eps)['T']-mkrow(i,k,2.0,PHI1,"phi1",eps)['lim']) for eps in EPS]
        out.append("| %s | %s, E=%g | %.3f | %.3f |"%("F (a)" if k==2 else "A (b)",dyn(i),S[i]['E'],math.log10(e[0]/e[1]),math.log10(e[1]/e[2])))
RATES="\n".join(out)
# P2: independence
out=[]
out.append("#### P2: phi- and c-independence, X_est = (T - V0)/(a - a') (F) resp. (T - V0)/(a'' - a) (A) -> X = C_1 (data)")
out.append("| (phi, c) | family | eps | state | T_eps | X_est | C_1 (data) | |X_est - C_1| |")
out.append("|---|---|---|---|---|---|---|---|")
cases=[("phi1, c=2","phi1",PHI1,2.0),("phi2, c=2","phi2",PHI2,2.0),("phi1, c=0.5","phi1",PHI1,0.5)]
for cn,phn,phi,cc in cases:
    for k in (2,1):
        for eps in (1e-3,1e-4):
            for i in (4,5,6):
                r=mkrow(i,k,cc,phi,phn,eps)
                if k==2: dd=a-a*(cc+phi[0])/(cc+phi[3])
                else: dd=a*(cc+phi[3])/(cc+phi[0])-a
                X=(r['T']-V0[i])/dd
                out.append("| %s | %s | %g | #%d | %s | %s | %s | %.1e |"%(cn,"F (a), k=2" if k==2 else "A (b), k=1",eps,i+1,f(r['T']),f(X),f(C1[i]),abs(X-C1[i])))
P2="\n".join(out)
# B1
out=[]
out.append("#### B1 control: A-site family at a = -1, c = -phi4 (kappa = z4): z1*T_eps = <B_1|Psi_eps> exactly")
out.append("| state | eps | z1 T_eps | <B_1|Psi_eps> | |difference| | C_1 (data) |")
out.append("|---|---|---|---|---|---|")
B1=prod_state([E(1),E(4)+E(1)]+[E(1)-E(4),E(4)+E(1)]*(Lp-1)).reshape(-1)
for i in (0,4,5,6):
    for eps in EPS:
        z=1+eps*PHI1
        r=mkrow(i,1,-PHI1[3],PHI1,"phi1",eps) if False else row(i,TW[("phi1",eps)][i][:4],1,-1.0,-PHI1[3],PHI1,eps,eps)
        d=B1@TW[("phi1",eps)][i][0]
        out.append("| %s, E=%g | %g | %s | %s | %.1e | %s |"%(dyn(i),S[i]['E'],eps,f(z[0]*r['T']),f(d),abs(z[0]*r['T']-d),f(C1[i])))
B1T="\n".join(out)
# P3: pair 2,3
out=[]
out.append("#### P3: mixed placements, phi1, c = 2, a = 1/3, eps = 1e-3 and 1e-4 (the same columns; k = 1, 2 repeated from P1 for the comparison at eps = 1e-4)")
out.append("| k (pair) | family | state | eps | T_eps | limit | |T - limit| | |LHS_prod - RHS_tau| |")
out.append("|---|---|---|---|---|---|---|---|")
for k in (1,2,3,4,5,6):
    for i in (0,4,5,6):
        for eps in (1e-3,1e-4):
            r=mkrow(i,k,2.0,PHI1,"phi1",eps)
            out.append("| %d (%d) | %s | %s, E=%g | %g | %s | %s | %.2e | %.1e |"%(k,(k+1)//2,"F (a)" if k%2==0 else "A (b)",dyn(i),S[i]['E'],eps,f(r['T']),f(r['lim']),abs(r['T']-r['lim']),abs(r['lhs_prod']-r['rhs_tau'])))
P3="\n".join(out)
pickle.dump(dict(P1=P1,RATES=RATES,P2=P2,B1=B1T,P3=P3),open('tables_part2.pkl','wb'))
for n,t in (("P1",P1),("RATES",RATES),("P2",P2),("B1",B1T),("P3",P3)): print(n,len(t))
