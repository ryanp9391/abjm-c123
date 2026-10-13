import numpy as np, time
from hw import *
Lp=3; N=6
THDIR=np.array(TH8[:N]); THABJM=np.array([1.,0,1,0,1,0])
PHI1=np.array([0.9,0.4,-0.3,-1.1]); PHI2=np.array([1.3,-0.2,0.5,-0.7])
SECS=sectors(Lp)
S=hwdata()['L3']
def prodv(sv): return prod_state(sv)
vac=[E(1),E(4)]*Lp
def bra_pair(p, A, F):  # vac everywhere except pair p (1-based) -> (A,F)
    sv=[E(1),E(4)]*Lp; sv[2*p-2]=A; sv[2*p-1]=F; return prod_state(sv)
V0=[prod_state(vac).reshape(-1)@s['psi'] for s in S]
X11=[[bra_pair(p,E(1),E(1)).reshape(-1)@s['psi'] for p in (1,2,3)] for s in S]
X44=[[bra_pair(p,E(4),E(4)).reshape(-1)@s['psi'] for p in (1,2,3)] for s in S]
C1=[s['C'][1] for s in S]
def twisted(idxs, th, z):
    """twisted ket eigenvectors for listed states (indices into S). returns dict i -> (psi_full, overlap)"""
    groups={}
    for i in idxs: groups.setdefault(tuple(S[i]['w']),[]).append(i)
    out={}
    for w,ii in groups.items():
        idx=SECS[w]; M=sector_matrix(Lp,idx,th,z)
        vals,vecs=np.linalg.eig(M); vecs=vecs/np.linalg.norm(vecs,axis=0)
        for i in ii:
            p0=S[i]['psi'][idx]
            ov=np.abs(vecs.conj().T@p0); k=np.argmax(ov); second=np.sort(ov)[-2] if len(ov)>1 else 0
            p=vecs[:,k]; dd=p0.conj()@p; p=p*abs(dd)/dd
            full=np.zeros(4**N,complex); full[idx]=p; out[i]=(full,ov[k],second,vals[k])
    return out
def dcl(k,th): return d_closed(k,th,Lp)
def ccl(k,th): return c_closed(k,th,Lp)
def site_vectors(k, a, th, z, eps, cc, phi, family=None):
    """product bra for the word whose last (free) letter is at site k.  sites j<k: torus rotated; site k: free family; sites >k: W_a
       returns sv list and the 'epsilon-normalisation' factor (eps*(c+phi_x))"""
    Lam=np.diag(z); LamI=np.diag(1/np.array(z)); I4=np.eye(4)
    sv=Wv(a)*Lp
    for j in range(k-1): sv=rot(sv,j,LamI if j%2==0 else Lam)
    if k%2==0:   # F site free: u (c eps - 1 + Lam)
        sv=rot(sv,k-1,(cc*eps-1)*I4+Lam); nrm=eps*(cc+phi[3])
    else:        # A site free: v (1 - kappa LamI)
        kap=1-cc*eps; sv=rot(sv,k-1,I4-kap*LamI); nrm=eps*(cc+phi[0])
    return sv,nrm
def limit_value(k, a, cc, phi, i, pair=None):
    p=(k+1)//2
    if k%2==0:
        ap=a*(cc+phi[0])/(cc+phi[3]); return V0[i]+(a-ap)*X11[i][p-1]
    else:
        app=a*(cc+phi[3])/(cc+phi[0]); return V0[i]+(app-a)*X44[i][p-1]
def tau(psi, th, z, u, Lbra=None):
    kk=np.argmax(np.abs(psi)); e=np.zeros(4**N); e[kk]=1
    y=t1_apply(e.reshape((4,)*N), th, z, u).reshape(-1)
    return (y@psi)/psi[kk]
def word_eig(psi, th, z, k, al, be):
    """eigenvalue prefactor of the word  prod_{j<k} t1(th_j)^{(-1)^j}  (al + be t1(th_k)^{+1 if k even else -1})"""
    ta=[tau(psi,th,z,th[j]) for j in range(k)]
    pre=1.0
    for j in range(1,k):  # letters 1..k-1: odd inverse, even forward
        pre*= (1/ta[j-1]) if j%2==1 else ta[j-1]
    last=ta[k-1]
    return pre*(al+be*(last if k%2==0 else 1/last))
def alpha_beta(k, a, cc, eps, th):
    if k%2==0: return (cc*eps-1)/ccl(k-1,th), 1/ccl(k,th)
    else:
        kap=1-cc*eps; return 1/ccl(k-1,th), -kap/ccl(k,th)
def row(i, tw, k, a, cc, phi, eps, delta):
    psi=tw[0]; th=THABJM+delta*THDIR; z=1+eps*phi
    sv,nrm=site_vectors(k,a,th,z,eps,cc,phi)
    bp=prod_state(sv).reshape(-1)
    Wb=prod_state(Wv(a)*Lp).reshape(-1)
    al,be=alpha_beta(k,a,cc,eps,th)
    ovW=Wb@psi
    lhs_prod=bp@psi
    rhs_tau=word_eig(psi,th,z,k,al,be)*ovW
    T=lhs_prod/nrm
    return dict(i=i,eps=eps,ov=tw[1],ovW=ovW,lhs_prod=lhs_prod,rhs_tau=rhs_tau,T=T,lim=limit_value(k,a,cc,phi,i),nrm=nrm)
if __name__=="__main__":
    pass
