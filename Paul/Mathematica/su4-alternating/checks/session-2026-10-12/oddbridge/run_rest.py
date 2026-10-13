import sys, time, pickle
sys.path.insert(0,'.')
from untw import *
np.set_printoptions(linewidth=200,precision=6)
a=1/3; cc=2.0; phi=PHI1
IDX=list(range(8,10))+list(range(10,11))+list(range(11,26))   # states 9..26 (0-based 8..25) -> skip big? include all
IDX=[i for i in range(26) if i not in (0,1,2,3,4,5,6)]
t0=time.time()
res={}
for eps in [1e-2,1e-3,1e-4]:
    th=THABJM+eps*THDIR; z=1+eps*phi
    tw=twisted(IDX,th,z)
    for i in IDX:
        for k in (1,2,3,4):
            r=row(i,tw[i],k,a,cc,phi,eps,eps)
            res[(i,k,eps)]=(abs(r['T']),abs(r['lim']),r['idres'] if 'idres' in r else abs(r['lhs_prod']-r['rhs_tau']),tw[i][1],tw[i][2])
    print("eps",eps,"%.1fs"%(time.time()-t0),flush=True)
pickle.dump(res,open('rest.pkl','wb'))
for i in IDX:
    print(i+1,S[i]['dyn'],"E=%.3f"%S[i]['E'],"|T| k=2: "+" ".join("%.2e"%res[(i,2,e)][0] for e in (1e-2,1e-3,1e-4)),"| k=1: "+" ".join("%.2e"%res[(i,1,e)][0] for e in (1e-2,1e-3,1e-4)),"idres max %.1e"%max(res[(i,k,e)][2] for k in (1,2,3,4) for e in (1e-2,1e-3,1e-4)),"ov %.4f 2nd %.4f"%(res[(i,2,1e-4)][3],res[(i,2,1e-4)][4]))
