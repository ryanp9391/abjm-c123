# Explicit ABJM 3-pt overlaps for on-shell highest-weight states, all bridge lengths, and the odd-bridge
# transfer-matrix formula through the untwisting limit.  Usage: python abjm_bridges.py [Lpairs]   (default 2)
# Conventions as README (A site rho(E_ij)=e_ij, F site rho(E_ij)=-e_ji, T = L_1...L_L, twist on the right).
# Bra of the note's 3-pt function, bridge l (fields): first l sites carry O2 fields (A: e1, F: e4),
# the rest O3 fields (A: e1-e4, F: e4+e1);  l = 0 is <W_{-1}^L|, l = 2L the vacuum.
import sys, itertools, functools, numpy as np
h=0.5; N=4
Lp=int(sys.argv[1]) if len(sys.argv)>1 else 2; L=2*Lp
def rA(i,j): m=np.zeros((N,N)); m[i,j]=1; return m
def rF(i,j): m=np.zeros((N,N)); m[j,i]=-1; return m
reps=[rA if s%2==0 else rF for s in range(L)]
I4=np.eye(N); e=[I4[i] for i in range(N)]
prod=lambda vs: functools.reduce(np.kron,vs)

def site_op(psi,s,M):                       # ket action of M on site s, psi of shape (4,)*L
    return np.moveaxis(np.tensordot(M,psi,axes=(1,s)),0,s)
def t1(u,psi,th,Lam):                       # t_1(u) psi = sum_ij T_ij(u) Lam_ji psi
    psi=psi.reshape((N,)*L).astype(complex); s=L-1
    Phi=[[h*site_op(psi,s,reps[s](k,j))+((u-th[s])*psi if k==j else 0) for j in range(N)] for k in range(N)]
    for s in range(L-2,-1,-1):
        Phi=[[ (u-th[s])*Phi[k][j]+h*sum(site_op(Phi[m][j],s,reps[s](k,m)) for m in range(N)) for j in range(N)] for k in range(N)]
    return sum(Phi[i][j]*Lam[j,i] for i in range(N) for j in range(N) if Lam[j,i]!=0).reshape(-1)
def Jtot(i,j,psi):
    psi=psi.reshape((N,)*L); return sum(site_op(psi,s,reps[s](i,j)) for s in range(L)).reshape(-1)
def shift(psi):                             # translation by one pair (2 sites)
    return np.moveaxis(np.moveaxis(psi.reshape((N,)*L),0,-1),0,-1).reshape(-1)

# weight sectors of product states
basis=list(itertools.product(range(N),repeat=L))
def weight(b):
    w=np.zeros(N,int)
    for s,a in enumerate(b): w[a]+= 1 if s%2==0 else -1
    return tuple(w)
sectors={}
for k,b in enumerate(basis): sectors.setdefault(weight(b),[]).append(k)
dominant=[w for w in sectors if all(w[i]>=w[i+1] for i in range(N-1))]
def unit(k): v=np.zeros(N**L); v[k]=1; return v

def bra(l,a=-1.0):
    O2=[e[0],e[3]]; O3=[e[0]+a*e[3], e[3]-a*e[0]]
    return prod([(O2 if s<l else O3)[s%2] for s in range(L)])

# ---- sanity: vector t1 vs matrix construction (L=4 sites only)
if L==4:
    import os; D=os.path.dirname(os.path.abspath(__file__))
    src=open(os.path.join(D,'g2test.py')).read().split('th=[0.2857')[0]; g={}; exec(src,g)
    thx=[2/7,-3/5,8/11,-11/8]; zx=np.array([1.5,1.25,2/3,4/7]); g['G']=np.diag(zx)
    Cm=g['Chain']([g['rA'],g['rF'],g['rA'],g['rF']],thx); psi=np.random.default_rng(0).standard_normal(256)
    print("vector t1 == matrix t1 (twisted, inhomogeneous):",np.allclose(Cm.t(1,0.37)@psi,t1(0.37,psi,thx,np.diag(zx))))

# ---- untwisted homogeneous chain: on-shell HW cyclic states
th0=np.zeros(L); Id=np.eye(N); us=(0.3711,-0.917)
states=[]
for w in sorted(dominant,reverse=True):
    idx=sectors[w]; B=np.array([unit(k) for k in idx]).T
    cons=[np.array([Jtot(i,i+1,B[:,c]) for c in range(B.shape[1])]).T for i in range(N-1)]
    cons.append(np.array([shift(B[:,c])-B[:,c] for c in range(B.shape[1])]).T)
    A=np.vstack(cons); sv=np.linalg.svd(A,compute_uv=False) if A.size else np.array([])
    _,s,Vh=np.linalg.svd(A); r=np.sum(s>1e-9); K=Vh[r:].conj().T
    if K.shape[1]==0: continue
    H=B@K
    M=np.array([[np.vdot(H[:,a],t1(us[0],H[:,b],th0,Id)+0.37*t1(us[1],H[:,b],th0,Id)) for b in range(H.shape[1])] for a in range(H.shape[1])])
    ev,V=np.linalg.eig(M)
    for k in range(len(ev)):
        psi=H@V[:,k]; psi=psi/np.linalg.norm(psi)
        tau=np.vdot(psi,t1(0.123,psi,th0,Id))
        gap=min([abs(ev[k]-ev[j]) for j in range(len(ev)) if j!=k],default=np.inf)
        states.append(dict(w=w,psi=psi,tau=tau,gap=gap))
dyn=lambda w: (w[0]-w[1],w[1]-w[2],w[2]-w[3])
print("\nL = %d pairs (%d sites): on-shell cyclic HW states: %d"%(Lp,L,len(states)))
print("Dynkin [p,q,r] | tau1(0.123) | min gap | C_l = <B_l|Psi>/|Psi|, l = 0..%d   (odd l = odd bridge)"%L)
for st in states:
    C=[bra(l)@st['psi'] for l in range(L+1)]
    st['C']=C
    print("%-10s %9.4f%+.4fi %8.1e | "%(str(list(dyn(st['w']))),st['tau'].real,st['tau'].imag,st['gap'])+"  ".join("%8.4f%+.4fi"%(c.real,c.imag) if abs(c.imag)>1e-9 else "%8.4f"%c.real for c in C))

# ---- odd-bridge formula through the untwisting limit (l = 1)
# <B_1|Psi> = (1 - z4/(c1 tau1(th1))) / (1 - z4/z1) * <W_{-1}^L|Psi>,  1/c1 = h prod_{F s} th1s prod_{A s>1}(th1s+h)
phi=np.array([0.9,0.4,-0.3,-1.1]); thdir=np.array([2/7,-3/5,8/11,-11/8,0.41,-0.77][:L])
def c1(th): return 1/(h*np.prod([th[0]-th[s] if s%2 else th[0]-th[s]+h for s in range(1,L)]))
print("\nOdd-bridge formula (l=1) on twisted, inhomogeneous deformations of each HW state:")
for st in states:
    if abs(st['C'][1])<1e-8: continue
    w=st['w']; idx=sectors[w]; B=np.array([unit(k) for k in idx]).T
    row=[]
    for eps,dl in ((1e-2,1e-2),(1e-3,1e-3),(1e-4,1e-4)):
        z=1+eps*phi; th=dl*thdir; Lam=np.diag(z)
        M=np.array([[B[:,a]@t1(us[0],B[:,b],th,Lam)+0.37*(B[:,a]@t1(us[1],B[:,b],th,Lam)) for b in range(B.shape[1])] for a in range(B.shape[1])])
        ev,V=np.linalg.eig(M); P=B@V; P=P/np.linalg.norm(P,axis=0)
        k=np.argmax(np.abs(P.conj().T@st['psi'])); psi=P[:,k]; psi=psi*np.vdot(psi,st['psi'])/abs(np.vdot(psi,st['psi']))
        tau=(t1(th[0],psi,th,Lam)@psi.conj())/(psi@psi.conj())
        lhs=bra(1)@psi; Wov=bra(0)@psi; fac=(1-z[3]/(c1(th)*tau))/(1-z[3]/z[0])
        row.append((eps,abs(lhs-fac*Wov)/abs(lhs),lhs,Wov,fac))
    print(" %-10s C_1(HW)=%s"%(str(list(dyn(w))),np.round(st['C'][1],6)))
    for eps,rel,lhs,Wov,fac in row:
        print("    eps=delta=%.0e: |LHS-RHS|/|LHS| = %.1e   LHS = %s   <W^L|Psi> = %.2e   bracket = %.2e"%(eps,rel,np.round(lhs,6),abs(Wov),abs(fac)))
