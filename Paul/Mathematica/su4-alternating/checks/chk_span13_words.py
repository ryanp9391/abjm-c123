import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'rotbase.py')).read().split('kap=0.37')[0])
vac=prod([e1,e4,e1,e4]); X,Y=T[0],T[2]
print("[X,Y]=0:",np.allclose(X@Y,Y@X))
mons=[(a,b) for a in range(3) for b in range(3)]
S=np.array([vac@np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Y,b) for a,b in mons]).T
sv=np.linalg.svd(S,compute_uv=False); print("dim span:",np.sum(sv>1e-10*sv[0]), " sv:",sv/sv[0])
# which site-2/4 content? project on V = span{e_i x e4 x e_j x e4}
Vb=np.array([prod([I[i],e4,I[j],e4]) for i in range(4) for j in range(4)]).T
Q,_=np.linalg.qr(S); Pv=Vb@Vb.T
M=Q.T@Pv@Q; ev=np.linalg.eigvalsh(M); print("dim(S cap V) =",np.sum(ev>1-1e-10), " eigen:",np.round(ev,8))
# basis of S cap V and the 4x4 matrices (site1 x site3) they give
w,U=np.linalg.eigh(M); K=Q@U[:,ev>1-1e-10]
mats=[ (Vb.T@K[:,k]).reshape(4,4) for k in range(K.shape[1])]
for m in mats: print(np.round(m/np.abs(m).max(),4))
print("S inside V (sites 2,4 = e4 x e4 for every word)?", np.linalg.norm(S-Pv@S)/np.linalg.norm(S))
# symmetric rank-1 targets v x v, v=e1 p(G^T), p quadratic: compare span of mats with {e1, e1G^T, e1G^T^2} sym products
B=[e1, e1@G.T, e1@G.T@G.T]; Ssym=np.array([(np.outer(B[i],B[j])+np.outer(B[j],B[i])).ravel()/2 for i in range(3) for j in range(i,3)]).T
Mm=np.array([m.ravel() for m in mats]).T
both=np.hstack([Mm,Ssym]); print("rank mats, rank Sym2(span e1 p(G^T)), rank joint:",np.linalg.matrix_rank(Mm,1e-9),np.linalg.matrix_rank(Ssym,1e-9),np.linalg.matrix_rank(both,1e-9))
def fit(tgt):
    c,*_=np.linalg.lstsq(S,tgt,rcond=None); return c, np.linalg.norm(S@c-tgt)/np.linalg.norm(tgt)
rng=np.random.default_rng(11)
print("words: ",["X^%dY^%d"%m for m in mons])
v=np.array([*rng.standard_normal(3),0]); c,r=fit(prod([v,e4,v,e4])); print("random v in span(e1,e2,e3), global v x e4 x v x e4: resid %.0e"%r)
k1,k2=0.37,-0.81; v=e1@(I-k1*G.T)@(I-k2*G.T); c,r=fit(prod([v,e4,v,e4])); print("v=e1(1-k1G^T)(1-k2G^T): resid %.0e  coeffs/c0:"%r, np.round(c/c[0],5))
v=e1@G.T@G.T; c,r=fit(prod([v,e4,v,e4])); print("pure (G^T)^2 on sites 1,3: resid %.0e"%r)
c,r=fit(prod([e1@G.T@G.T,e4,e1,e4])); print("(G^T)^2 on site 1 alone: resid %.0e"%r)
v=np.array([0.3,1,-0.7,0.4]); c,r=fit(prod([v,e4,v,e4])); print("v with v4!=0 (not global): resid %.0e"%r)
