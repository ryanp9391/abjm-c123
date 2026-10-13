import sys, time, pickle
sys.path.insert(0,'.')
from untw import *
np.set_printoptions(linewidth=200,precision=6)
IDX=[0,1,2,3,4,5,6,7,8,9,10]   # states 1-11 (sector <= 42)
EPS=[1e-2,1e-3,1e-4]
t0=time.time()
TW={}
for phn,phi in (("phi1",PHI1),("phi2",PHI2)):
    for eps in EPS:
        th=THABJM+eps*THDIR; z=1+eps*phi
        TW[(phn,eps)]=twisted(IDX,th,z)
        print("twisted",phn,eps,"%.1fs"%(time.time()-t0),flush=True)
pickle.dump(TW,open('tw_cache.pkl','wb'))
