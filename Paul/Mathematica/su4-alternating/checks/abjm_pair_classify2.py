import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'abjm_pair_classify.py')).read().split('s1=["1"]')[0])
import time; t0=time.time()
s1=["1"]+[n for n in names if "th1" in n]; s2=["1"]+[n for n in names if "th2" in n]
P1={}; 
for i,a_ in enumerate(s1):
    for b_ in s1[i:]:
        P1[a_+"."+b_]=W@op(a_)@op(b_)
P2={}
for i,a_ in enumerate(s2):
    for b_ in s2[i:]:
        P2[a_+"."+b_]=op(a_)@op(b_)
print(len(P1),len(P2))
found={}
for na,x in P1.items():
    if not np.isfinite(x).all() or np.linalg.norm(x)<1e-12: continue
    for nb,Mb in P2.items():
        y=x@Mb
        if not np.isfinite(y).all(): continue
        ny=np.linalg.norm(y)
        if ny<1e-12: continue
        y=y/ny
        if V12(y)>1e-9: continue
        M=pairstate(y)
        if isinstance(M,str): continue
        e=str(expo(M)); found.setdefault(e,[]).append(na+" | "+nb)
for k,v in sorted(found.items(),key=lambda kv:len(kv[1])): print(k,": %d words, e.g."%len(v),v[:4])
print("time %.0fs"%(time.time()-t0))
