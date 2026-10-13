import os, tempfile
import itertools, numpy as np, functools, collections
from fractions import Fraction
# variable layout: v2:0-3, u2:4-7, v3:8-11, u3:12-15
def comps(n,k=4):
    if k==1: yield (n,); return
    for i in range(n+1):
        for r in comps(n-i,k-1): yield (i,)+r
def weight(e):
    w=[0]*4
    for i in range(4):
        w[i]=e[i]-e[4+i]+e[8+i]-e[12+i]
    return tuple(w)
def monos(n2,n3,w,bps):
    out=[]
    for a in comps(n2):
        for b in comps(n2):
            for c in comps(n3):
                for d in comps(n3):
                    e=a+b+c+d
                    if weight(e)!=tuple(w): continue
                    if bps and ((e[0]>0 and e[4]>0) or (e[8]>0 and e[12]>0)): continue
                    out.append(e)
    return out
@functools.lru_cache(maxsize=None)
def reduce_mono(e):
    """normal form of monomial modulo v2.u2, v3.u3 (leading v_{a1}u_{a1} -> -sum_{i>=2} v_i u_i). returns dict"""
    e=list(e)
    for (iv,iu) in ((0,4),(8,12)):
        if e[iv]>0 and e[iu]>0:
            res=collections.defaultdict(int)
            base=e[:]; base[iv]-=1; base[iu]-=1
            for i in range(1,4):
                f=base[:]; f[iv+i]+=1; f[iu+i]+=1
                for k,c in reduce_mono(tuple(f)).items(): res[k]-=c
            return {k:c for k,c in res.items() if c}
    return {tuple(e):1}
def Dop(e,i,j,bps=True):
    """D_ij on monomial e: sum_v v_i d/dv_j - sum_u u_j d/du_i ; returns dict"""
    res=collections.defaultdict(int)
    for off in (0,8):   # vectors v2,v3
        k=off+j
        if e[k]>0:
            f=list(e); f[k]-=1; f[off+i]+=1
            res[tuple(f)]+=e[k]
    for off in (4,12): # covectors u2,u3
        k=off+i
        if e[k]>0:
            f=list(e); f[k]-=1; f[off+j]+=1
            res[tuple(f)]-=e[k]
    out=collections.defaultdict(int)
    for f,c in res.items():
        if bps:
            for g,cc in reduce_mono(f).items(): out[g]+=c*cc
        else: out[f]+=c
    return {k:c for k,c in out.items() if c}
def inv_dim(n2,n3,w,bps):
    basis=monos(n2,n3,w,bps)
    if not basis: return 0,0
    rows=[]
    blocks=[]
    for i in range(3):
        tgt={}
        cols=[]
        for e in basis:
            d=Dop(e,i,i+1,bps)
            cols.append(d)
            for k in d: tgt.setdefault(k,len(tgt))
        A=np.zeros((len(tgt),len(basis)))
        for c,d in enumerate(cols):
            for k,v in d.items(): A[tgt[k],c]=v
        blocks.append(A)
    A=np.vstack(blocks)
    rank=np.linalg.matrix_rank(A) if A.size else 0
    return len(basis)-rank,len(basis)
if __name__=='__main__':
    import pickle
    from core import S
    from inv import monomials
    res=pickle.load(open(os.path.join(tempfile.gettempdir(),'twoblock_fits.pkl'),'rb'))
    import sys
    Ls=[int(x) for x in sys.argv[1:]] or [2,3]
    seen=set()
    bad=[]
    for L in Ls:
        for k,s in enumerate(S[L]):
            for M in range(L+1):
                for bps in (True,False):
                    key=(L,M,tuple(int(x) for x in s['w']),bps)
                    if key in seen: continue
                    seen.add(key)
                    dim,nb=inv_dim(L-M,M,s['w'],bps)
                    nm=res[(L,k,M,bps)]['nm']
                    if dim!=nm: bad.append((key,dim,nm))
    print('sectors tested',len(seen),'mismatches',bad)
