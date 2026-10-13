import json, numpy as np, sympy as sp, time
from core import S, DIG
d=json.load(open(r'C:\Users\pryan\git-projects\abjm-c123\Paul\Mathematica\su4-alternating\data\abjm-twoblock-forms.json'))
v2=sp.symbols('v2_1:5'); u2=sp.symbols('u2_1:5'); v3=sp.symbols('v3_1:5'); u3=sp.symbols('u3_1:5')
allv=list(v2)+list(u2)+list(v3)+list(u3)
I=dict(a2=v2[0],a3=v3[0],d2=u2[3],d3=u3[3],mv=v2[0]*v3[1]-v2[1]*v3[0],mu=u2[2]*u3[3]-u2[3]*u3[2],
       p22=sum(a*b for a,b in zip(v2,u2)),p23=sum(a*b for a,b in zip(v2,u3)),p32=sum(a*b for a,b in zip(v3,u2)),p33=sum(a*b for a,b in zip(v3,u3)))
NAMES=['a2','a3','d2','d3','mv','mu','p22','p23','p32','p33']
def Fsym(L,M,psi):
    D=DIG[L]; n=2*L
    nz=np.nonzero(np.abs(psi)>0)[0]
    tot=0
    for Iidx in nz:
        t=sp.Float(psi[Iidx].real)+sp.I*sp.Float(psi[Iidx].imag)
        for s in range(n):
            k=s//2; blk=0 if k<L-M else 1; typ=s%2
            vec=[[v2,u2],[v3,u3]][blk][typ]
            t*=vec[D[Iidx,s]]
        tot+=t
    return sp.expand(tot)
def Dder(F,i,j):
    r=0
    for v in (v2,v3): r+=v[i]*sp.diff(F,v[j])
    for u in (u2,u3): r-=u[j]*sp.diff(F,u[i])
    return sp.expand(r)
t0=time.time()
worst={2:[0,0],3:[0,0]}
for L in (2,3):
    for k,s in enumerate(S[L]):
        for M in range(L+1):
            F=Fsym(L,M,s['psi'])
            if F==0: continue
            P=sp.Poly(F,*allv)
            # derivation check
            for i in range(3):
                dd=sp.Poly(Dder(F,i,i+1),*allv)
                if dd.terms(): worst[L][0]=max(worst[L][0],max(abs(complex(c)) for _,c in dd.terms()))
            G=0
            for t in d['L%d'%L][k]['M'][str(M)]['full']['terms']:
                c=complex(t['c'][0],t['c'][1])
                m=1
                for n_,e in zip(NAMES,t['exp']): m*=I[n_]**e
                G+=(sp.Float(c.real)+sp.I*sp.Float(c.imag))*m
            diff=sp.Poly(sp.expand(F-G),*allv)
            if diff.terms(): worst[L][1]=max(worst[L][1],max(abs(complex(c)) for _,c in diff.terms()))
    print('sympy L',L,'max |D_{i,i+1}F| coefficient',worst[L][0],'max |F - sum c I| coefficient',worst[L][1],time.time()-t0)
