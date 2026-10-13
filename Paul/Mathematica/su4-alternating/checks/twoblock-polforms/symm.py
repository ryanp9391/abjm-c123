import os, tempfile
import pickle, numpy as np, collections
from core import S
res=pickle.load(open(os.path.join(tempfile.gettempdir(),'twoblock_fits.pkl'),'rb'))
TH=1e-9
def coeffs(L,k,M,bps=True):
    r=res[(L,k,M,bps)]
    d={}
    if r['nm']==0: return d
    for m,c in zip(r['mons'],r['c']):
        d[tuple(m.get(n,0) for n in ['a2','a3','d2','d3','mv','mu','p22','p23','p32','p33'])]=c
    return d
def swap(key):
    a2,a3,d2,d3,mv,mu,p22,p23,p32,p33=key
    return (a3,a2,d3,d2,mv,mu,p33,p32,p23,p22)
def sgn(key): return (-1)**(key[4]+key[5])
# 1. cyclicity: F_M(2,3) = F_{L-M}(3,2): c(M;key)=sgn(key)*c(L-M;swap(key))
for bps in (True,False):
  worst=0;cnt=0
  for L in (2,3,4):
    for k in range(len(S[L])):
        for M in range(L+1):
            c1=coeffs(L,k,M,bps); c2=coeffs(L,k,L-M,bps)
            keys=set(c1)|set(swap(x) for x in c2)
            for key in keys:
                a=c1.get(key,0); b=sgn(key)*c2.get(swap(key),0)
                worst=max(worst,abs(a-b)); cnt+=1
  print('bps',bps,'cyclicity (M <-> L-M with 2<->3 swap, mv,mu -> -mv,-mu): worst diff',worst,'over',cnt)
# 2. note's C: F at note polarisation
e=np.eye(4)
v2,u2,v3,u3=e[0],e[3],e[0]-e[3],e[3]+e[0]
a2,a3,d2,d3=v2[0],v3[0],u2[3],u3[3]
mv=v2[0]*v3[1]-v2[1]*v3[0]; mu=u2[2]*u3[3]-u2[3]*u3[2]
p22,p23,p32,p33=v2@u2,v2@u3,v3@u2,v3@u3
print('note pol: a2,a3,d2,d3,mv,mu,p22,p23,p32,p33 =',a2,a3,d2,d3,mv,mu,p22,p23,p32,p33)
vals=dict(zip(['a2','a3','d2','d3','mv','mu','p22','p23','p32','p33'],[a2,a3,d2,d3,mv,mu,p22,p23,p32,p33]))
worst=0
for L in (2,3,4):
    for k,s in enumerate(S[L]):
        for M in range(L+1):
            c=coeffs(L,k,M,True)
            val=0
            for key,ci in c.items():
                val+=ci*np.prod([vals[n]**ee for n,ee in zip(['a2','a3','d2','d3','mv','mu','p22','p23','p32','p33'],key)])
            worst=max(worst,abs(val-s['C'][2*(L-M)]))
print('C_{2(L-M)} (saved) vs invariant form at note polarisation: worst diff',worst)
