# Negative check: naive anti-fundamental-auxiliary Lax u - (+/-) h sum e_ij x rho(E_ji) does NOT commute with t1, t3.
import numpy as np, os, functools
os.chdir(os.path.dirname(os.path.abspath(__file__)))
g={}; exec(open('g2test.py').read().split('th=[0.2857')[0],g)
I=np.eye(4); h=0.5; reps=[g['rA'],g['rF'],g['rA'],g['rF']]
def so(a,M): return functools.reduce(np.kron,[M if b==a else I for b in range(4)])
def tbar(u,th,sgn):
    T=[[np.eye(256)*(i==j) for j in range(4)] for i in range(4)]
    for a,r in enumerate(reps):
        La=[[so(a,(u-th[a])*(i==j)*I+sgn*h*r(j,i)) for j in range(4)] for i in range(4)]
        T=[[sum(T[i][k]@La[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    return sum(T[i][i] for i in range(4))
for th in ([0,0,0,0],[2/7,-3/5,8/11,-11/8]):
    g['G']=np.eye(4); C=g['Chain'](reps,th); t1=C.t(1,0.3711); t3=C.t(3,0.41)
    for sgn in (-1,1):
        tb=tbar(0.29,th,sgn)
        print("th=%s sgn=%+d: [t1,tbar]=0 %s  [t3,tbar]=0 %s"%(th[:2],sgn,np.allclose(t1@tb,tb@t1),np.allclose(t3@tb,tb@t3)))
