import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'scan_blocks.py')).read().split('k1,k2=0.37')[0])
def classify(refs,psi,Om):
    if np.linalg.norm(psi)<1e-10*np.linalg.norm(Om): return "zero"
    c=np.vdot(Om,psi)/np.vdot(Om,Om)
    if np.linalg.norm(psi-c*Om)<1e-9*np.linalg.norm(psi): return "prop"
    T4=psi.reshape(4,4,4,4); ranks=[]
    for s in range(4):
        sv=np.linalg.svd(np.moveaxis(T4,s,0).reshape(4,-1),compute_uv=False); ranks.append(int(np.sum(sv>1e-9*sv[0])))
    return "product" if max(ranks)==1 else "entangled ranks %s"%ranks
def scan2(label,refs,want):
    Om=prod(refs); vecs={n:Om@ops[n] for n in names}; rows={}
    words=[(n,vecs[n]) for n in names]+[(a+" . "+b,vecs[a]@ops[b]) for i,a in enumerate(names) for b in names[i:]]
    for w,psi in words:
        cl=classify(refs,psi,Om)
        if cl in("zero","prop"): continue
        for B in loc_blocks(refs,psi):
            if B in want: rows.setdefault(B,[]).append((w,cl))
    print("=====",label)
    for B in want:
        L=rows.get(B,[]); print(" sites",[b+1 for b in B],": %d words"%len(L))
        singles=[x for x in L if " . " not in x[0]]; prs=[x for x in L if " . " in x[0]]
        for w,cl in singles: print("     %-22s %s"%(w,cl))
        for w,cl in prs[:6]: print("     %-40s %s"%(w,cl))
        if len(prs)>6: print("     ... (+%d pairs; product ones: %d)"%(len(prs)-6,sum(c=='product' for _,c in prs)))
want=[(0,),(1,),(2,),(3,),(0,1),(2,3),(1,2),(0,2),(1,3),(0,3),(0,1,2),(1,2,3)]
k1,k2=0.37,-0.81; v=e1@(I-k1*G.T)@(I-k2*G.T)
scan2("Omega_2",[v,e4,v,e4],want)
rng=np.random.default_rng(1); gr=np.eye(4)+0.4*rng.standard_normal((4,4))
scan2("generic global rotation",[e1@gr,e4@np.linalg.inv(gr).T]*2,want)
