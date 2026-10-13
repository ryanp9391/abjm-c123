import numpy as np, io, contextlib, collections
with contextlib.redirect_stdout(io.StringIO()):
    from inv import *
Pstar=[np.array(x,dtype=complex) for x in ([1,0,0,0],[0,0,0,1],[1,1,0,1],[1,1,0,-2])]
Pnote=[np.array(x,dtype=complex) for x in ([1,0,0,0],[0,0,0,1],[1,0,0,-1],[1,0,0,1])]
for P in (Pstar,Pnote):
    I=invs(*P); print({k:(complex(v)) for k,v in I.items()})
    x=I['p23']/(I['a2']*I['d3']); y=I['p32']/(I['a3']*I['d2']); muv=I['mv']/(I['a2']*I['a3']); muu=I['mu']/(I['d2']*I['d3'])
    print('x,y,muv,muu',x,y,muv,muu)
vis=collections.Counter(); tot=collections.Counter()
for L in (2,3,4):
    for k,s in enumerate(S[L]):
        for M in range(1,L):
            Fn=Fval(s['psi'],L,M,Pnote)/ (1.0)
            Fs=Fval(s['psi'],L,M,Pstar)
            Fvs=(1)**0
            # Fvac at P*: a2^(L-M) d2^(L-M) a3^M d3^M
            I=invs(*Pstar); Fv=(I['a2']*I['d2'])**(L-M)*(I['a3']*I['d3'])**M
            tot[L]+=1
            if abs(Fn)<1e-12 and abs(Fs)>1e-9:
                vis[(L,tuple(int(z) for z in s['dyn']),M)]+=1
print(vis)
