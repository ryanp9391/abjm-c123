import sys; sys.path.insert(0,'.')
from untw import *
np.set_printoptions(linewidth=200,precision=6)
# forward word identities at twisted ABJM-point parameters (L=3, eps=delta=1e-3) and at the standard point
def word_check(k, a, cc, eps, th, z, phi):
    sv,nrm=site_vectors(k,a,th,z,eps,cc,phi)
    bp=prod_state(sv)
    Wb=prod_state(Wv(a)*Lp)
    al,be=alpha_beta(k,a,cc,eps,th)
    T=lambda x,j: t1_apply(x,th,z,th[j-1])
    if k%2==0:
        lhs=bp
        for j in range(1,k,2): lhs=T(lhs,j)
        rhs=Wb
        for j in range(2,k-1,2): rhs=T(rhs,j)
        rhs=al*rhs+be*T(rhs,k)
    else:
        lhs=bp
        for j in range(1,k+1,2): lhs=T(lhs,j)
        rhs=Wb
        for j in range(2,k,2): rhs=T(rhs,j)
        rhs=al*T(rhs,k)+be*rhs
    return np.linalg.norm(lhs-rhs)/np.linalg.norm(rhs)
for eps in [1e-3]:
    th=THABJM+eps*THDIR; z=1+eps*PHI1
    print("ABJM-point twisted eps=%g:"%eps,[ "%.1e"%word_check(k,1/3,2.0,eps,th,z,PHI1) for k in range(1,7)])
th=np.array(TH8[:6]); z=np.array(Z4)
print("standard point (a=1/3, c=2, eps=0.05 with phi1):",[ "%.1e"%word_check(k,1/3,2.0,0.05,th,np.array(Z4),PHI1) for k in range(1,7)])
# B1 identity at eps
eps=1e-3; th=THABJM+eps*THDIR; z=1+eps*PHI1
sv,nrm=site_vectors(1,-1.0,th,z,eps,-PHI1[3],PHI1)
bp=prod_state(sv).reshape(-1)
B1=prod_state([E(1),E(4)+E(1)]+[E(1)-E(4),E(4)+E(1)]*(Lp-1)).reshape(-1)
print("bra_prod vs (1-z4/z1) B1:", np.linalg.norm(bp-(1-z[3]/z[0])*B1))
# random BPS pair: chain constants at L=3
rng=np.random.default_rng(1)
v=rng.integers(-5,6,4).astype(float); r=rng.integers(-5,6,4).astype(float); u=r-(r@v)/(v@v)*v
print("random BPS v.u=",v@u)
th=np.array(TH8[:6]); z=np.array(Z4)
sv0=[v,u]*3; Lam=np.diag(z); LamI=np.diag(1/z)
def P(k):
    sv=list(sv0)
    for s in range(k): sv=rot(sv,s,LamI if s%2==0 else Lam)
    return prod_state(sv)
errs=[]
for k in range(1,7):
    if k%2==1: y=t1_apply(P(k),th,z,th[k-1]); r_,dev,zc=ratio(y,P(k-1)); dk=1/r_
    else: y=t1_apply(P(k-1),th,z,th[k-1]); r_,dev,zc=ratio(y,P(k)); dk=r_
    errs.append((abs(dk-dcl(k,th))/abs(dcl(k,th)),dev))
print("random BPS pair L=3: (relerr d_k, prod dev):",["%.1e/%.1e"%e for e in errs])
