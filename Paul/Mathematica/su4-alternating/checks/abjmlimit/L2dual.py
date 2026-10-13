import sympy as sp, pickle
from itertools import permutations, combinations
from chain import perm_sign
exec(open('L2num2.py').read().split('results = []')[0])
eps = sp.symbols('eps')
keys = [(a,k) for a in range(4) for k in range(4)]
# generic-d words M(eps) along d = 2h + eps
def Mof(dval):
    Y = tnum(1, 0, dval); rows=[]
    for a in range(4):
        v = vac*tnum(a, dval-h, dval)
        for k in range(4): rows.append(v); v = v*Y
    return sp.Matrix.vstack(*rows)
Me = Mof(2*h+eps).applyfunc(sp.expand)
Minv = Me.inv()  # rational in eps
reg = pickle.load(open('L2reg0.pkl','rb')); R = reg['R']; Ri = R.inv()
def is_reg(a,k): return a==3 or (a==2 and k<=2)
# Laurent of columns of Minv
out = []
for i,x in enumerate(keys):
    col = Minv[:, i]
    lead = col.applyfunc(lambda e: sp.limit(eps*sp.cancel(e), eps, 0))
    fin = None
    order = 1 if any(e != 0 for e in lead) else 0
    if order:
        # residue equals which dual-of-R column?
        match = [keys[j] for j in range(16) if (lead - Ri[:, j]).is_zero_matrix]
        matchneg = [keys[j] for j in range(16) if (lead + Ri[:, j]).is_zero_matrix]
        out.append((x, 'pole 1/eps', 'residue = +r_%s' % (match,) if match else ('residue = -r_%s' % (matchneg,) if matchneg else 'other')))
    else:
        lim = col.applyfunc(lambda e: sp.limit(sp.cancel(e), eps, 0))
        match = [keys[j] for j in range(16) if (lim - Ri[:, j]).is_zero_matrix]
        out.append((x, 'finite', 'limit = r_%s' % (match,) if match else 'limit != a single dual'))
for o in out: print(o)
# finite parts of the (1,k+1) columns: Minv[:, (1,k+1)] + Minv[:, (2,k)] -> ?
for k in range(3):
    s = (Minv[:, keys.index((1,k+1))] + Minv[:, keys.index((2,k))]).applyfunc(lambda e: sp.limit(sp.cancel(e), eps, 0))
    print('lim [y_(1,%d) + y_(2,%d)] = r_(1,%d)?' % (k+1, k, k+1), (s - Ri[:, keys.index((1,k+1))]).is_zero_matrix)
pickle.dump({'Ri': Ri}, open('L2dual.pkl','wb'))
