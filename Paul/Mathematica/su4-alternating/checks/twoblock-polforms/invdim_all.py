import itertools
from invdim import *
from inv import monomials
tot=0;bad=[]
for L in (2,3,4):
    for M in range(L+1):
        n2,n3=L-M,M
        for bps in (True,False):
            for w in itertools.product(range(-L,L+1),repeat=3):
                w=w+(-sum(w),)
                if abs(w[3])>L: continue
                dim,nb=inv_dim(n2,n3,w,bps)
                _,mons=monomials(L,M,w,bps)
                tot+=1
                if dim!=len(mons): bad.append((L,M,bps,w,dim,len(mons)))
print('all weights tested',tot,'bad',bad[:10],len(bad))
