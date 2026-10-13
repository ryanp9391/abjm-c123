import numpy as np, itertools
from core import *
rng=np.random.default_rng(1)
NAMES=['a2','a3','d2','d3','mv','mu','p22','p23','p32','p33']
def invs(v2,u2,v3,u3):
    return dict(a2=v2[0],a3=v3[0],d2=u2[3],d3=u3[3],
        mv=v2[0]*v3[1]-v2[1]*v3[0], mu=u2[2]*u3[3]-u2[3]*u3[2],
        p22=v2@u2,p23=v2@u3,p32=v3@u2,p33=v3@u3)
def randpt(bps=True):
    v2,u2,v3,u3=[rng.normal(size=4)+1j*rng.normal(size=4) for _ in range(4)]
    v2,u2,v3,u3=[x/np.linalg.norm(x) for x in (v2,u2,v3,u3)]
    if bps:
        for v,u in ((v2,u2),(v3,u3)):
            w=np.conj(v)
            u-= (v@u)/(v@w)*w
        assert abs(v2@u2)<1e-12 and abs(v3@u3)<1e-12
    return v2,u2,v3,u3
def Fval(psi,L,M,pt):
    v2,u2,v3,u3=pt
    vecs=[]
    for k in range(L):
        if k<L-M: vecs+= [v2,u2]
        else: vecs+=[v3,u3]
    return overlap(psi,L,vecs)
def monomials(L,M,w,bps):
    n2,n3=L-M,M
    out=[]
    names=NAMES if not bps else [x for x in NAMES if x not in('p22','p33')]
    # exponents
    e5=w[1]; e6=-w[2]
    if e5<0 or e6<0: return names,[]
    rng_=range(0,L+1)
    for e1 in range(0,L+1):
      for e2 in range(0,L+1):
        if e1+e2+e5!=w[0]: continue
        for e7 in (range(0,L+1) if not bps else [0]):
         for e8 in range(0,L+1):
          # v2: e1+e5+e7+e8=n2
          if e1+e5+e7+e8!=n2: continue
          for e9 in range(0,L+1):
           for e10 in (range(0,L+1) if not bps else [0]):
            if e2+e5+e9+e10!=n3: continue
            # u2: e3+e6+e7+e9=n2 ; u3: e4+e6+e8+e10=n3
            e3=n2-e6-e7-e9; e4=n3-e6-e8-e10
            if e3<0 or e4<0: continue
            if -(e3+e4+e6)!=w[3]: continue
            out.append(dict(a2=e1,a3=e2,d2=e3,d3=e4,mv=e5,mu=e6,p22=e7,p23=e8,p32=e9,p33=e10))
    return names,out
def fit(L,s,M,bps=True,npts=None):
    names,mons=monomials(L,M,s['w'],bps)
    nm=len(mons)
    npts=npts or max(40,3*nm+10)
    pts=[randpt(bps) for _ in range(npts)]
    F=np.array([Fval(s['psi'],L,M,p) for p in pts])
    if nm==0:
        return dict(nm=0,rank=0,res=np.abs(F).max()/ (np.abs(F).max()+1e-300), Fmax=np.abs(F).max())
    A=np.empty((npts,nm),dtype=complex)
    for r,p in enumerate(pts):
        I=invs(*p)
        for c,m in enumerate(mons):
            A[r,c]=np.prod([I[k]**m[k] for k in names if m.get(k,0)])
    sv=np.linalg.svd(A,compute_uv=False)
    rank=int((sv>sv[0]*1e-9).sum())
    c,*_=np.linalg.lstsq(A,F,rcond=1e-10)
    res=np.abs(A@c-F).max()
    return dict(nm=nm,rank=rank,res=res,Fmax=np.abs(F).max(),c=c,mons=mons,names=names)
if __name__=='__main__':
    for L in (2,3):
        for k,s in enumerate(S[L]):
            for M in range(L+1):
                r=fit(L,s,M,True)
                print(L,k,s['dyn'],M,'nmon',r['nm'],'rank',r['rank'],'Fmax %.3g res %.2g'%(r['Fmax'],r['res']))
