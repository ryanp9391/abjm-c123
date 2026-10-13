import sympy as sp, mpmath as mp, pickle, time
mp.mp.dps = 60
d, uu = sp.symbols('d uu')
h = sp.Rational(1,2)
tsym = pickle.load(open('L2tsym.pkl','rb'))
tsym[0] = sp.eye(16)
vac = sp.zeros(1,16); vac[0,3] = 1
def tnum(a, uval, dval):
    return tsym[a].subs({uu: uval, d: dval}) if a else sp.eye(16)
def words(dval):
    Y = tnum(1, 0, dval); out = {}
    for a in range(4):
        v = vac*tnum(a, dval-h, dval)
        for k in range(4):
            out[(a,k)] = v; v = v*Y
    return out
keys = [(a,k) for a in range(4) for k in range(4)]
# test bra: product v (A) x w (F)
vA = sp.Matrix([[1, sp.Rational(2,3), sp.Rational(-1,2), sp.Rational(3,5)]])
wF = sp.Matrix([[sp.Rational(2,7), -1, sp.Rational(1,3), 1]])
b = sp.kronecker_product(vA, wF)
def tomp(Mx): return mp.matrix([[mp.mpf(sp.Rational(x).p)/sp.Rational(x).q for x in Mx.row(i)] for i in range(Mx.rows)])
def eigstates(dval):
    O = tomp(tnum(1, sp.Rational(3137,10000), dval))
    E, ER = mp.eig(O)
    return [ER[:, i] for i in range(16)]
def tau(a, uval, dval, psi):
    Mt = tomp(tnum(a, uval, dval)) if a else mp.eye(16)
    v = Mt*psi
    j = max(range(16), key=lambda i: abs(psi[i]))
    return v[j]/psi[j]
def dtau(a, uval, dval, psi):
    Mt = tomp(sp.diff(tsym[a], uu).subs({uu: uval, d: dval}))
    v = Mt*psi
    j = max(range(16), key=lambda i: abs(psi[i]))
    return v[j]/psi[j]
results = []
for delta in [sp.Rational(1,100), sp.Rational(1,1000), sp.Rational(1,10000)]:
    dval = 2*h + delta
    W = words(dval)
    M = sp.Matrix.vstack(*[W[k] for k in keys])
    c = (b * M.inv())  # exact
    psis = eigstates(dval)
    maxc = max(abs(float(x)) for x in c)
    res = []
    for psi in psis[:16]:
        taus = {a: tau(a, dval-h, dval, psi) for a in range(4)}
        tF = tau(1, 0, dval, psi)
        S = mp.fsum(mp.mpf(sp.Rational(c[i]).p)/sp.Rational(c[i]).q * taus[k[0]]*tF**k[1] for i,k in enumerate(keys))
        bpsi = mp.fsum(mp.mpf(sp.Rational(b[i]).p)/sp.Rational(b[i]).q*psi[i] for i in range(16))
        direct = bpsi/psi[3]
        res.append(abs(S-direct)/max(1,abs(direct)))
    # residues
    resid = {k: float(delta*c[i]) for i,k in enumerate(keys) if k[0]==3 or (k[0]==2 and k[1]<=2)}
    sums = {k: float(c[keys.index((2,k))] + c[keys.index((1,k+1))]) for k in range(3)}
    print('delta', delta, 'max|c|', maxc, ' max rel |sum c tau - <b|psi>/<0|psi>| over 16 states', mp.nstr(max(res),3), flush=True)
    print('   delta*c (a=2,k<=2; a=3):', {k: round(v,8) for k,v in resid.items()})
    print('   c(2,k)+c(1,k+1):', {k: round(v,8) for k,v in sums.items()})
pickle.dump(None, open('L2num2.done','wb'))
