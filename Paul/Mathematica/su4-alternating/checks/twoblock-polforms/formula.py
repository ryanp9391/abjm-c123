import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from symm import res
from core import S
def nstruct(L,M,p,q,r):
    n2,n3=L-M,M
    if (2*q+r-p)%4 or (p+2*q-r)%4: return None
    e5=(2*q+r-p)//4; e6=(p+2*q-r)//4
    if e5<0 or e6<0: return 0
    s2=(p+2*q+r)
    if s2%2: return None
    s=L-s2//2
    if s<0: return 0
    lo=max(0,s-(n3-e5),s-(n2-e6)); hi=min(s,n2-e5,n3-e6)
    return max(0,hi-lo+1)
bad=0;tot=0
for L in (2,3,4):
    for k,s in enumerate(S[L]):
        p,q,r=[int(x) for x in s['dyn']]
        for M in range(L+1):
            tot+=1
            if nstruct(L,M,p,q,r)!=res[(L,k,M,True)]['nm']: bad+=1; print('mismatch',L,k,(p,q,r),M,nstruct(L,M,p,q,r),res[(L,k,M,True)]['nm'])
print('formula checked for',tot,'(state,M); mismatches',bad)
