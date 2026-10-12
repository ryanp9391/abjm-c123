# L = 2: resolve the degenerate [1,0,1] HW cyclic pair of abjm_bridges.py with the fused t3 and print its C_l.
import numpy as np, sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv=['x','2']; _f=__file__; __file__=os.path.abspath('abjm_bridges.py')
exec(open('abjm_bridges.py').read().split('# ---- odd-bridge formula')[0])
g={}; exec(open('g2test.py').read().split('th=[0.2857')[0],g)
g['G']=np.eye(4); Cm=g['Chain']([g['rA'],g['rF'],g['rA'],g['rF']],[0,0,0,0])
t3=Cm.t(3,0.41)+0.3*Cm.t(3,-0.77); t1m=Cm.t(1,0.3711)
print("[t1,t3]=0:",np.allclose(t1m@t3,t3@t1m))
deg=[s for s in states if [int(x) for x in dyn(s['w'])]==[1,0,1]]
Hs=np.array([s['psi'] for s in deg]).T
M3=np.linalg.pinv(Hs)@t3@Hs; print("t3 on the degenerate [1,0,1] pair:\n",np.round(M3,6))
ev,V=np.linalg.eig(M3); print("t3 eigenvalues:",ev)
for k in range(2):
    psi=Hs@V[:,k]; psi=psi/np.linalg.norm(psi)
    print("  state %d: C_l ="%k,[complex(np.round(bra(l)@psi,6)) for l in range(5)])
