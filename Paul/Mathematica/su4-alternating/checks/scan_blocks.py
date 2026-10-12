import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'rotbase.py')).read().split('kap=0.37')[0])
import pickle, itertools as it, collections
import tempfile; ops=pickle.load(open(os.path.join(tempfile.gettempdir(),'su4_alt_ops.pkl'),'rb')); names=list(ops)  # run scan_ops_build.py first
blocks=[b for r in (1,2,3) for b in it.combinations(range(4),r)]
def loc_blocks(refs,psi):
    T4=psi.reshape(4,4,4,4); nrm=np.linalg.norm(psi); out=[]
    for B in blocks:
        x=T4
        for s in range(4):
            if s in B: continue
            r=refs[s]/np.linalg.norm(refs[s]); P=np.eye(4)-np.outer(r,r)   # component orthogonal to ref on site s
            x2=np.moveaxis(np.tensordot(P,np.moveaxis(T4,s,0),axes=(1,0)),0,s)
            if np.linalg.norm(x2)>1e-9*nrm: break
        else: out.append(B)
    # minimal blocks only
    return [B for B in out if not any(set(C)<set(B) for C in out)]
def scan(label,refs,pairs=True):
    Om=prod(refs); res=collections.defaultdict(list)
    vecs={n:Om@ops[n] for n in names}
    for n in names:
        for B in loc_blocks(refs,vecs[n]): res[B].append(n)
    if pairs:
        for i,n1 in enumerate(names):
            for n2 in names[i:]:
                for B in loc_blocks(refs,vecs[n1]@ops[n2]): res[B].append(n1+" . "+n2)
    print("=====",label)
    for B in blocks:
        L=res.get(B,[])
        if L:
            singles=[x for x in L if " . " not in x]; prs=[x for x in L if " . " in x]
            print(" sites",[b+1 for b in B],": %d singles, %d pairs"%(len(singles),len(prs)))
            print("     singles:",singles[:14])
            if prs: print("     pairs e.g.:",prs[:8])
k1,k2=0.37,-0.81
v=e1@(I-k1*G.T)@(I-k2*G.T)
scan("bare <0|",[e1,e4,e1,e4])
scan("Omega_2 = v x e4 x v x e4, v=e1(1-k1G^T)(1-k2G^T)",[v,e4,v,e4])
