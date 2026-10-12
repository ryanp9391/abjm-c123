import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'abjm_pair_base.py')).read().split('def target')[0])
import itertools as it
def tgt(ap,rest=None): # pair1 = W_{a'}, pair2 = W_a
    x,y=pairW(ap); return prod([x,y,vA,wF])
basis3=[prod([e1,e4,vA,wF]), prod([e4,e4,vA,wF])-prod([e1,e1,vA,wF]), prod([e4,e1,vA,wF])]
def inspan(S,t): c,*_=np.linalg.lstsq(S,t,rcond=None); return np.linalg.norm(S@c-t)/np.linalg.norm(t)
print("16 SoV words: W_{a'} family basis resid:",["%.0e"%inspan(S,b) for b in basis3])
# enlarge: library of t_k(th_s + m h) for s=1,2 and products of two
lib={}
for s in (0,1):
    for m in (-2,-1,0,1,2):
        for k in (1,2,3):
            M=C4.t(k,th[s]+m*h); lib["t%d(th%d%+dh)"%(k,s+1,m)]=M
            if np.linalg.cond(M)<1e12: lib["t%d(th%d%+dh)^-1"%(k,s+1,m)]=np.linalg.inv(M)
names=list(lib); v1={n:W@lib[n] for n in names}
zero=[n for n in names if np.linalg.norm(v1[n])<1e-12]; print('ops killing W:',zero)
def clean(cs):
    cs=[c/np.linalg.norm(c) for c in cs if np.isfinite(c).all() and np.linalg.norm(c)>1e-12]
    return cs
cols=clean([W]+list(v1.values())); S1=np.array(cols).T
print("singles(sites1,2 pts): rank",np.linalg.matrix_rank(S1,1e-9),"; family basis resid:",["%.0e"%inspan(S1,b) for b in basis3])
cols2=cols+[v1[a]@lib[b] for i,a in enumerate(names) for b in names[i:]]
S2=np.array(clean(cols2)).T
print("pairs: rank",np.linalg.matrix_rank(S2,1e-9),"; family basis resid:",["%.0e"%inspan(S2,b) for b in basis3])
