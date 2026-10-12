import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'rotbase.py')).read().split('kap=0.37')[0])
z=np.array([1.5,1.25,2/3,4/7]); Lam=np.diag(z); G=Lam      # diagonal (physical) frame twist
C4=Chain([rA,rF,rA,rF],th)
pairW=lambda a:(e1+a*e4, e4-a*e1)
def Wstate(avals): return prod([x for a in avals for x in pairW(a)])
a=-1.0; vA,wF=pairW(a); W=Wstate([a,a])
V12=lambda psi: np.linalg.norm(psi.reshape(16,16)-np.outer(psi.reshape(16,16)@np.kron(vA,wF)/np.dot(np.kron(vA,wF),np.kron(vA,wF)),np.kron(vA,wF)))/np.linalg.norm(psi)
tA={k:C4.t(k,th[0]-h) for k in range(4)}; tF=C4.t(1,th[1])
words={}
for k in range(4):
    for n in range(4): words["t%d(th1-h)t1(th2)^%d"%(k,n)]=tA[k]@np.linalg.matrix_power(tF,n)
vecs={n:W@M for n,M in words.items()}
print("words localising on pair 1 (sites 3,4 stay W):",[n for n,v in vecs.items() if V12(v)<1e-10])
S=np.array(list(vecs.values())).T
# subspace of span(S) lying in V12: solve S c = x (x) Wpair  ->  project out components outside V12
pw=np.kron(vA,wF); pw=pw/np.linalg.norm(pw); Pout=np.kron(np.eye(16),np.eye(16)-np.outer(pw,pw))
ns=np.linalg.svd(Pout@S)[2]; sv=np.linalg.svd(Pout@S,compute_uv=False)
null=ns[np.sum(sv>1e-9*sv[0]):]; print("dim span(words) =",np.linalg.matrix_rank(S,1e-9*np.linalg.norm(S)), "  dim(span cap V12) =",null.shape[0])
def target(kap):
    D=np.diag(1-kap*z); return prod([vA@D, wF@np.linalg.inv(D), vA, wF])
for kap in (0.3,-1.1):
    t=target(kap); c,*_=np.linalg.lstsq(S,t,rcond=None); r=np.linalg.norm(S@c-t)/np.linalg.norm(t)
    print("kappa=%.4f  rho_12(1-kappa Lam) W in span: resid %.1e"%(kap,r))
