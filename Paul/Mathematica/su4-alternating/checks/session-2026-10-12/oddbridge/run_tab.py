import sys, time, pickle
sys.path.insert(0,'.')
from untw import *
np.set_printoptions(linewidth=200,precision=6)
TW=pickle.load(open('tw_cache.pkl','rb'))
EPS=[1e-2,1e-3,1e-4]
def fmt(x):
    x=complex(x); return "%.5g%+.5gj"%(x.real,x.imag) if abs(x.imag)>1e-13 else "%.6g"%x.real
a=1/3
cases=[("phi1,c=2",PHI1,"phi1",2.0),("phi2,c=2",PHI2,"phi2",2.0),("phi1,c=0.5",PHI1,"phi1",0.5)]
for k in range(1,7):
    fam="F-site (a)" if k%2==0 else "A-site (b)"
    print("=== position k=%d (pair %d) %s"%(k,(k+1)//2,fam))
    for cn,phi,phn,cc in cases:
        print("  case",cn)
        for i in [0,4,5,6,1,2,7]:
            line=[]
            for eps in EPS:
                r=row(i,TW[(phn,eps)][i],k,a,cc,phi,eps,eps)
                line.append((eps,abs(r['T']-r['lim']),abs(r['lhs_prod']-r['rhs_tau']),r['T'],r['lim'],r['ov']))
            print("   state %2d %s: "%(i+1,S[i]['dyn'])+" | ".join("eps=%g err=%.2e idres=%.1e"%(e,er,ir) for e,er,ir,_,_,_ in line)+"  T(1e-4)=%s lim=%s"%(fmt(line[-1][3]),fmt(line[-1][4])))
