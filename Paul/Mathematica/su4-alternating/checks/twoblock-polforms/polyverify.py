import numpy as np, collections, io, contextlib, pickle, itertools, time
with contextlib.redirect_stdout(io.StringIO()):
    from symm import coeffs, res
from core import S, DIG
from invdim import reduce_mono, Dop, weight
# polynomial as dict {exp16: complex}
def Fpoly(psi,L,M):
    D=DIG[L]; n=2*L
    nz=np.nonzero(np.abs(psi)>0)[0]
    out=collections.defaultdict(complex)
    for I in nz:
        e=[0]*16
        for s in range(n):
            k=s//2   # pair index
            blk=0 if k<L-M else 2   # block 2 -> vars 0..7 ; block 3 -> 8..15
            typ=s%2                 # 0: A (v) ; 1: F (u)
            base=blk*4+typ*4 if blk==0 else 8+typ*4
            e[base+D[I,s]]+=1
        out[tuple(e)]+=psi[I]
    return out
def pmul(p,q):
    out=collections.defaultdict(complex)
    for a,ca in p.items():
        for b,cb in q.items():
            out[tuple(x+y for x,y in zip(a,b))]+=ca*cb
    return out
def ppow(p,n):
    r={tuple([0]*16):1.0}
    for _ in range(n): r=pmul(r,p)
    return r
def var(i):
    e=[0]*16; e[i]=1; return {tuple(e):1.0}
def padd(p,q,cq=1.0):
    out=collections.defaultdict(complex,p)
    for k,v in q.items(): out[k]+=cq*v
    return out
def inv_polys():
    v2=[var(i) for i in range(0,4)]; u2=[var(i) for i in range(4,8)]
    v3=[var(i) for i in range(8,12)]; u3=[var(i) for i in range(12,16)]
    def dot(a,b):
        r={}
        for x,y in zip(a,b): r=padd(r,pmul(x,y))
        return r
    I={}
    I['a2']=v2[0]; I['a3']=v3[0]; I['d2']=u2[3]; I['d3']=u3[3]
    I['mv']=padd(pmul(v2[0],v3[1]),pmul(v2[1],v3[0]),-1.0)
    I['mu']=padd(pmul(u2[2],u3[3]),pmul(u2[3],u3[2]),-1.0)
    I['p22']=dot(v2,u2); I['p23']=dot(v2,u3); I['p32']=dot(v3,u2); I['p33']=dot(v3,u3)
    return I
INV=inv_polys()
NAMES=['a2','a3','d2','d3','mv','mu','p22','p23','p32','p33']
_cache={}
def mono_poly(key):
    if key in _cache: return _cache[key]
    r={tuple([0]*16):1.0}
    for n,e in zip(NAMES,key):
        if e: r=pmul(r,ppow(INV[n],e))
    _cache[key]=r
    return r
def Dpoly(p,i,j):
    out=collections.defaultdict(complex)
    for e,c in p.items():
        for f,cc in Dop(e,i,j,bps=False).items():
            out[f]+=c*cc
    return out
def maxabs(p): return max([abs(v) for v in p.values()]+[0])
def reduce_poly(p):
    out=collections.defaultdict(complex)
    for e,c in p.items():
        for f,cc in reduce_mono(e).items(): out[f]+=c*cc
    return out
if __name__=='__main__':
    import sys
    Ls=[int(x) for x in sys.argv[1:]] or [2,3]
    for L in Ls:
        worstD=0; worstW=0; worstE=0; worstB=0; nterms=[]; t0=time.time()
        for k,s in enumerate(S[L]):
            for M in range(L+1):
                F=Fpoly(s['psi'],L,M)
                F={e:c for e,c in F.items() if abs(c)>0}
                # Borel invariance for all i<j
                for i in range(4):
                    for j in range(i+1,4):
                        worstD=max(worstD,maxabs(Dpoly(F,i,j)))
                # torus weight
                wset=set(weight(e) for e,c in F.items() if abs(c)>1e-14)
                if wset and wset!={tuple(int(x) for x in s['w'])}: worstW=max(worstW,1)
                # expansion
                r=res[(L,k,M,False)]
                G={}
                if r['nm']:
                    for m,c in zip(r['mons'],r['c']):
                        key=tuple(m.get(n,0) for n in NAMES)
                        G=padd(G,mono_poly(key),c)
                diff=padd(F,G,-1.0)
                worstE=max(worstE,maxabs(diff))
                # BPS: reduce both mod ideal
                rb=res[(L,k,M,True)]
                Gb={}
                if rb['nm']:
                    for m,c in zip(rb['mons'],rb['c']):
                        key=tuple(m.get(n,0) for n in NAMES)
                        Gb=padd(Gb,mono_poly(key),c)
                dB=padd(reduce_poly(F),reduce_poly(Gb),-1.0)
                worstB=max(worstB,maxabs(dB))
                nterms.append((len([1 for c in F.values() if abs(c)>1e-14]),len([1 for c in reduce_poly(F).values() if abs(c)>1e-14])))
        print('L',L,'max|D_ij F|',worstD,'weight mismatch',worstW,'max|F - sum c I| (full 16-var polynomial)',worstE,
              'max |red(F) - red(sum_BPS c I)|',worstB,'time',time.time()-t0)
        print('   max terms in F (no BPS) / normal form (BPS):',max(a for a,b in nterms),max(b for a,b in nterms))
