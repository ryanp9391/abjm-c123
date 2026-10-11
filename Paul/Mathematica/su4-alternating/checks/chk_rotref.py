# Globally rotated reference states v(x)w(x)v(x)w: localisation criterion w G v^T = 0, iterated common rotations.
import os,sys
D=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(D,'chk7.py')).read().split('r,c=prop(vac@Y')[0]
src=src.replace("exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'g2test.py')).read()","exec(open(os.path.join(D,'g2test.py')).read()")
exec(src)
np.set_printoptions(precision=4,suppress=True,linewidth=150)
T=[C4.t(1,th[s]) for s in range(4)]; Ti=[np.linalg.inv(x) for x in T]
prod=lambda vs: functools.reduce(np.kron,vs)
def local(state_vecs,out):
    # is out in span{ reference with site k replaced by anything }, for some k? return per-site residual
    res=[]
    for k in range(4):
        B=np.array([prod([I[j] if b==k else state_vecs[b] for b in range(4)]) for j in range(4)]).T
        c,*_=np.linalg.lstsq(B,out,rcond=None); res.append(np.linalg.norm(B@c-out)/np.linalg.norm(out))
    return res
def report(name,vecs):
    print("==",name, " v.w =",[round(float(vecs[0]@vecs[1]),12)])
    Om=prod(vecs)
    for nm,op in (("t1(th1)^-1",Ti[0]),("t1(th2)",T[1]),("t1(th3)",T[2]),("t1(th4)^-1",Ti[3]),("t1(th1)",T[0]),("t1(th4)",T[3])):
        r=local(vecs,Om@op); print("  %-11s per-site resid"%nm, ["%.0e"%x for x in r])
    # {1,3} common rotation g=1-kappa G^T on A sites ; {2,4} common rotation on F sites with plain G (site vec w -> w(1-kappa G))
    for lab,(i,j),M in (("{1,3} g=1-kG^T",(0,2),G.T),("{2,4} 1-kG",(1,3),G)):
        ws=[Om,Om@T[i],Om@T[j],Om@T[i]@T[j]]
        worst=0
        for kap in (0.3,-1.7):
            nv=list(vecs); nv[i]=vecs[i]@(I-kap*M); nv[j]=vecs[j]@(I-kap*M); tgt=prod(nv)
            A=np.array(ws).T; c,*_=np.linalg.lstsq(A,tgt,rcond=None); worst=max(worst,np.linalg.norm(A@c-tgt)/np.linalg.norm(tgt))
        print("  %-16s in span{1,t_i,t_j,t_it_j}: resid %.0e"%(lab,worst))
kap=0.37; c1=G[0,0]
gp=I-kap*(c1*e(0,0)+e(0,1))
print("F vac fixed by g':",np.allclose(e4@np.linalg.inv(gp).T,e4), "  e1 g' == e1 (1-kG^T):",np.allclose(e1@gp,e1@(I-kap*G.T)))
report("bare <0|",[e1,e4,e1,e4])
v=e1@gp; report("Omega = <0|rho(g'), g'=exp E12-type",[v,e4,v,e4])
gq=I-0.53*(e(2,3)+0.4*e(3,3)) # F-side analogue: acts as gq^{-T} on F
wq=e4@np.linalg.inv(gq).T; print("A vac fixed by gq:",np.allclose(e1@gq,e1), " wq=",wq)
report("<0|rho(gq), E34-type",[e1,wq,e1,wq])
report("<0|rho(g' gq)",[e1@gp@gq,e4@np.linalg.inv(gp@gq).T]*2)
v2=v@(I-0.8*G.T); report("second {1,3} step: v(1-kG^T), still global? v4=%.2g"%v2[3],[v2,e4,v2,e4])
v3=v2@(I-0.6*G.T); report("third step",[v3,e4,v3,e4])
rng=np.random.default_rng(1); gr=np.eye(4)+0.4*rng.standard_normal((4,4))
report("generic global rotation",[e1@gr,e4@np.linalg.inv(gr).T]*2)
# ---- iterated common rotations; criterion on random good references
def Pcoef(vecs,i,j,M,kap):
    Om=prod(vecs); ws=[Om@T[i],Om@T[j],Om@T[i]@T[j]]
    nv=list(vecs); nv[i]=vecs[i]@(I-kap*M); nv[j]=vecs[j]@(I-kap*M); tgt=prod(nv)
    A=np.array(ws+[-tgt]).T; s,*_=np.linalg.lstsq(A,-Om,rcond=None)
    a,b,c,lam=s; P=np.eye(256)+a*T[i]+b*T[j]+c*T[i]@T[j]
    return P/lam, np.linalg.norm(Om@P/lam-tgt)/np.linalg.norm(tgt), (a/lam,b/lam,c/lam,1/lam)
bare=[e1,e4,e1,e4]; vac=prod(bare)
k1,k2=0.37,-0.81
P13,r,cf=Pcoef(bare,0,2,G.T,k1); print("P13 on <0|: resid %.0e coeffs(1,a,b,c)*norm"%r,np.array(cf))
v=e1@(I-k1*G.T); Om=[v,e4,v,e4]
_,r,cfO=Pcoef(Om,0,2,G.T,k1); print("P13 on Omega (same kappa): coeffs",np.array(cfO))
# two common {1,3} steps
P13b,_,_=Pcoef(bare,0,2,G.T,k2)
s=vac@P13@P13b; vv=e1@(I-k1*G.T)@(I-k2*G.T); tgt=prod([vv,e4,vv,e4])
print("<0|P13(k1)P13(k2) = (e1(1-k1G^T)(1-k2G^T) x e4)^2 ? resid %.0e   v=%s"%(np.linalg.norm(s-tgt)/np.linalg.norm(tgt),vv))
# {1,3} x {2,4}: one global rotation on all four sites
P24,r,cf24=Pcoef(bare,1,3,G,k2); print("P24 on <0|: resid %.0e"%r)
s=vac@P13@P24; v=e1@(I-k1*G.T); w=e4@(I-k2*G); tgt=prod([v,w,v,w])
print("<0|P13 P24 = (v x w)^2: resid %.0e ; v.w=%.1e ; v=%s w=%s"%(np.linalg.norm(s-tgt)/np.linalg.norm(tgt),v@w,v,w))
# global rotation n with v=row1, w=row4 of n^{-T}: find n
# criterion test: random v, w with v.w=0 and w G v^T=0
rng=np.random.default_rng(7)
for trial in range(2):
    v=rng.standard_normal(4); Mx=np.array([v,G@v]); ns=np.linalg.svd(Mx)[2][2:]; w=rng.standard_normal(2)@ns
    report("random good ref (v.w=0, wGv=0) #%d"%trial,[v,w,v,w])
v=rng.standard_normal(4); ns=np.linalg.svd(v[None])[2][1:]; w=rng.standard_normal(3)@ns
print("wGv for generic:",w@G@v); report("random global ref with wGv!=0",[v,w,v,w])
