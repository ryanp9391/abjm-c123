import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'abjm_pair_classify2.py')).read().split('found={}')[0])
locw={}
for na,x in P1.items():
    if not np.isfinite(x).all() or np.linalg.norm(x)<1e-12: continue
    for nb,Mb in P2.items():
        y=x@Mb
        if not np.isfinite(y).all() or np.linalg.norm(y)<1e-12: continue
        y=y/np.linalg.norm(y)
        if V12(y)<1e-9: locw[na+" | "+nb]=y
L=np.array(list(locw.values())).T
print("words localising on pair 1:",len(locw)," span dim:",np.linalg.matrix_rank(L,1e-9))
S12=prod([np.zeros(4),np.zeros(4),vA,wF])
sing=(np.kron(np.eye(4).reshape(-1),np.kron(vA,wF))); sing/=np.linalg.norm(sing)
for nm,t in (("vacuum pair",basis3[0]),("family b2",basis3[1]),("family b3",basis3[2]),("singlet_12 x W",sing),("site-2 alone Lam",prod([vA,wF@Lam,vA,wF]))):
    c,*_=np.linalg.lstsq(L,t,rcond=None); print("  %-18s in span: resid %.0e"%(nm,np.linalg.norm(L@c-t)/np.linalg.norm(t)))
# what does the localised span look like on pair 1?  4x4 matrices M (site1 x site2)
B=np.linalg.svd(L,full_matrices=False)[0][:,:np.linalg.matrix_rank(L,1e-9)]
print("pair-1 support of localised span (rows=site1 comp, cols=site2 comp), max |entry| over basis:")
mx=np.zeros((4,4))
for k in range(B.shape[1]): mx=np.maximum(mx,np.abs((B[:,k].reshape(16,16)@np.kron(vA,wF)).reshape(4,4)))
print(np.round(mx,3))
