import sympy as sp, time, pickle
from L2sym import *
eps = sp.symbols('eps')
# (a) d-line: d = 2h + eps, u_A = theta_A - h
th_e = [2*h+eps, 0]
Ae = [t(a, 2*h+eps-h, th_e) for a in range(4)]
Ye = t(1, 0, th_e)
# (a') fixed theta, u-shift
th0 = [2*h, 0]
Ap = [t(a, 2*h-h+eps, th0) for a in range(4)]
Y0 = t(1, 0, th0)
def words2(A, Y):
    out = {}
    for a in range(4):
        v = vac*A[a]
        for k in range(4):
            out[(a,k)] = v.applyfunc(sp.expand); v = v*Y
    return out
We = words2(Ae, Ye); Wp = words2(Ap, Y0)
def regbasis(W):
    reg = {}
    for (a,k), v in W.items():
        if a == 3:
            r = (v/eps)
        elif a == 2 and k <= 2:
            r = (v - W[(1,k+1)])/eps
        else:
            r = v
        reg[(a,k)] = r.applyfunc(lambda x: sp.limit(sp.cancel(x), eps, 0) if x != 0 else 0)
    return reg
t0=time.time()
Re = regbasis(We); Rp = regbasis(Wp)
keys = sorted(Re.keys())
MRe = sp.Matrix.vstack(*[Re[k] for k in keys]); MRp = sp.Matrix.vstack(*[Rp[k] for k in keys])
print('rank reg (a) d-derivative:', MRe.rank(), ' rank reg (a\') u-derivative:', MRp.rank(), time.time()-t0)
# compare: same vectors for unregularised ones obviously; regularised ones differ?
for k in keys:
    if (k[0]==3) or (k[0]==2 and k[1]<=2):
        diff = (Re[k]-Rp[k]).applyfunc(sp.nsimplify)
        # is difference in span of the 9 unregularised words?
        base = sp.Matrix.vstack(*[Re[q] for q in keys if not ((q[0]==3) or (q[0]==2 and q[1]<=2))])
        inspan = sp.Matrix.vstack(base, diff).rank() == base.rank()
        ratio = None
        print(k, 'd-deriv - u-deriv in span of order-0 words:', inspan, ' d-deriv proportional to u-deriv:', sp.Matrix.vstack(Re[k],Rp[k]).rank()==1)
pickle.dump({'Re':Re,'Rp':Rp,'We':We,'Wp':Wp}, open('L2reg.pkl','wb'))
