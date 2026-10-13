import sympy as sp, mpmath as mp, pickle
exec(open('L2num2.py').read().split('results = []')[0])
d0 = 2*h; u0 = d0 - h
Y = tnum(1, 0, d0)
T1p = sp.diff(tsym[1], uu).subs({uu: u0, d: d0}); T2p = sp.diff(tsym[2], uu).subs({uu: u0, d: d0}); T3p = sp.diff(tsym[3], uu).subs({uu: u0, d: d0})
Wop = T2p - T1p*Y; Zop = T3p
A0 = [tnum(a, u0, d0) for a in range(4)]
print('t3(u0)=0:', A0[3] == sp.zeros(16,16), ' t2(u0)=t1(u0)t1(thF):', A0[2] == A0[1]*Y, ' [W, t1(u)] = 0:', (Wop*tnum(1, sp.Rational(1,7), d0) - tnum(1, sp.Rational(1,7), d0)*Wop) == sp.zeros(16,16))
def regop(a,k):
    if a < 2 or (a == 2 and k == 3): return A0[a]
    return Wop if a == 2 else Zop
R = sp.Matrix.vstack(*[vac*regop(a,k)*Y**k for (a,k) in keys])
print('rank R =', R.rank())
r = b*R.inv()
print('regularised coefficients r (a=2,k<=2; a=3):', {k: float(r[i]) for i,k in enumerate(keys) if k[0]==3 or (k[0]==2 and k[1]<=2)})
print('r(1,k+1) for k=0..2:', [float(r[keys.index((1,k+1))]) for k in range(3)])
psis = eigstates(d0)
errs = []
for psi in psis:
    t1 = tau(1, u0, d0, psi); t2 = tau(2, u0, d0, psi); tF = tau(1, 0, d0, psi)
    w = dtau(2, u0, d0, psi) - dtau(1, u0, d0, psi)*tF; y = dtau(3, u0, d0, psi)
    lam = {0: 1, 1: t1, 2: w, 3: y}
    S = mp.fsum(mp.mpf(sp.Rational(r[i]).p)/sp.Rational(r[i]).q * (t2 if (k[0]==2 and k[1]==3) else lam[k[0]]) * tF**k[1] for i,k in enumerate(keys))
    bpsi = mp.fsum(mp.mpf(sp.Rational(b[i]).p)/sp.Rational(b[i]).q*psi[i] for i in range(16))
    errs.append(abs(S - bpsi/psi[3])/max(1, abs(bpsi/psi[3])))
print('d=2h: max rel |sum r_j lambda_j - <b|psi>/<0|psi>| over 16 eigenstates:', mp.nstr(max(errs), 3))
pickle.dump({'r': r, 'R': R}, open('L2reg0.pkl','wb'))
