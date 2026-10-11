import os
import numpy as np, functools, itertools
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ffbar2.py')).read().split('# L=2: diagonal')[0])
th=[0.2857,-0.6,0.731,-1.37]; reps=[rf,rb,rf,rb]
t1,so=chain(reps,th)
I=np.eye(N)
Kop=lambda a,b: sum(so(a,e(i,j))@so(b,e(i,j)) for i in range(N) for j in range(N))
Pop=lambda a,b: sum(so(a,e(i,j))@so(b,e(j,i)) for i in range(N) for j in range(N))
Ks={(a,b):Kop(a,b) for a in range(4) for b in range(4) if a!=b}
Ps={(a,b):Pop(a,b) for a in range(4) for b in range(4) if a!=b}
Id=np.eye(N**4)
Gv={'G':G,'GT':G.T,'Gi':np.linalg.inv(G),'GiT':np.linalg.inv(G).T}
for a in range(4):
    M=t1(th[a]); others=[b for b in range(4) if b!=a]; found=False
    for types in itertools.product(['+K','-K','+P','-P'],repeat=3):
        Rs={}
        for b,tp in zip(others,types):
            X=Ks[(a,b)] if tp[1]=='K' else Ps[(a,b)]
            Rs[b]=(th[a]-th[b])*Id+(1 if tp[0]=='+' else -1)*h*X
        for perm in itertools.permutations(others):
            for gpos in range(4):
                for gn,g in Gv.items():
                    ops=[Rs[b] for b in perm]; ops.insert(gpos,so(a,g))
                    Y=functools.reduce(np.matmul,ops)
                    c=np.vdot(Y.ravel(),M.ravel())/np.vdot(Y.ravel(),Y.ravel())
                    if np.linalg.norm(c*Y-M)/np.linalg.norm(M)<1e-9:
                        print("site",a+1,"c=%.4g"%c.real,"R types",dict(zip([o+1 for o in others],types)),"order",[ ('G'+gn if x is None else x+1) for x in (list(perm[:gpos])+[None]+list(perm[gpos:]))]); found=True
                        break
                if found: break
            if found: break
        if found: break
    if not found: print("site",a+1,"no qKZ fit")
