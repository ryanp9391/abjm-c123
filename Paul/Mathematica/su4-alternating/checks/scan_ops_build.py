import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'rotbase.py')).read().split('kap=0.37')[0])
import time, pickle, itertools as it
t0=time.time()
ops={}
for s in range(4):
    for m in (-2,-1,0,1,2):
        u=th[s]+m*h
        for a in (1,2,3):
            M=C4.t(a,u); ops[("t%d(th%d%+d h)"%(a,s+1,m))]=M
            if np.linalg.cond(M)<1e12: ops[("t%d(th%d%+d h)^-1"%(a,s+1,m))]=np.linalg.inv(M)
print("ops:",len(ops)," time %.0fs"%(time.time()-t0))
import tempfile; pickle.dump(ops,open(os.path.join(tempfile.gettempdir(),'su4_alt_ops.pkl'),'wb'))
