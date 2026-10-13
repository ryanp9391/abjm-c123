import numpy as np, sys, time
from ob import *
import t_const
from t_const import run
rng=np.random.default_rng(11)
# L=4 at standard point
t0=time.time()
res=run(4, TH8, Z4)
print("L=4 standard (%.1fs)"%(time.time()-t0))
for k,dk,dc,dev,zc in res: print("  k=%d d=%.12g closed=%.12g relerr=%.2e prodDev=%.2e zeroComp=%.2e"%(k,dk,dc,abs(dk-dc)/abs(dc),dev,zc))
# random points
for Lp in [2,3]:
    for trial in range(2):
        th=list(rng.uniform(-1.5,1.5,2*Lp)); z=list(rng.uniform(0.4,2.0,4)); a=rng.uniform(-1,1)
        res=run(Lp, th, z, a)
        print("L=%d random point: max relerr %.2e, max prodDev %.2e"%(Lp,max(abs(r[1]-r[2])/abs(r[2]) for r in res),max(r[3] for r in res)))
# cumulative c at L=2
th=TH8[:4]
print("c3,c4 at L=2:", c_closed(3,th,2), c_closed(4,th,2), "expected -321399232/112262625 =", -321399232/112262625, " -1108723/174000=", -1108723/174000)
