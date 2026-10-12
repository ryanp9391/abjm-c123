import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'abjm_pair_span.py')).read().split('cols2=cols')[0])
import itertools as it
op=lambda n: np.eye(256) if n=="1" else lib[n]
# exponent fit: pair1 = (d1 e1 + a d4 e4) x (f4 e4 - a f1 e1), pair2 = W ; identify d,f as monomials in z
def pairstate(x):
    if V12(x)>1e-9: return None
    M=(x.reshape(16,16)@np.kron(vA,wF)).reshape(4,4)/np.dot(np.kron(vA,wF),np.kron(vA,wF))
    if np.linalg.matrix_rank(M,1e-9*np.abs(M).max())!=1: return "entangled"
    return M
def expo(M):
    # M = u v^T with u=(d1,0,0,a d4), v=(-a f1,0,0,f4); return monomial exponents (i,j) with d=z^i (A), f=z^j (F) if they fit
    c=[M[0,3],M[0,0],M[3,3],M[3,0]]   # d1f4, -a d1 f1, a d4 f4, -a^2 d4 f1
    v=np.array([c[0],-c[1]/a,c[2]/a,-c[3]/a**2]); v=v/v[1]   # normalise by d1 f1
    for i in range(-3,4):
        for j in range(-3,4):
            w=np.array([z[0]**i*z[3]**j, z[0]**i*z[0]**j, z[3]**i*z[3]**j, z[3]**i*z[0]**j]); w=w/w[1]
            if np.allclose(v,w,rtol=1e-8): return (i,j)
    return ("other",np.round(v,6))
s1=["1"]+[n for n in names if "th1" in n]; s2=["1"]+[n for n in names if "th2" in n]
found={}
for a_ in s1:
    x=W@op(a_)
    for b_ in s2:
        y=x@op(b_)
        if not np.isfinite(y).all() or np.linalg.norm(y)<1e-12: continue
        M=pairstate(y/np.linalg.norm(y))
        if M is None or isinstance(M,str): continue
        e=expo(M); found.setdefault(str(e),[]).append(a_+" * "+b_)
for k,v in found.items(): print(k,":",v[:6],"(+%d)"%(len(v)-6) if len(v)>6 else "")
