# torus rotations rho(X) = X (x) X^{-T}, X = f(G^T): equals a polynomial in the SoV words as an OPERATOR
import sympy as sp, pickle
from chain import *
U = sp.Symbol('U')
hN = sp.Rational(1, 2); th1N, th2N = sp.Rational(2, 7), sp.Rational(-3, 5)
zN = [sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(4, 7)]
G = companion(esym(zN))
tops = pickle.load(open('step2.pkl', 'rb'))['tops']
thS, hS = sp.symbols('theta h')
cc = pickle.load(open('step1c.pkl', 'rb'))['coef']
vS = sp.symbols('v1:5'); uS = sp.symbols('u1:5')
tA = {0: sp.eye(16)}; tA.update({a: tops[a].subs(U, th1N - hN) for a in (1, 2, 3)})
t1F = tops[1].subs(U, th2N)
words = {(a, k): tA[a] * t1F**k for a in range(4) for k in range(4)}
# commutation check
print('t_a(th1-h), t1(th2) commute:', all((tA[a] * t1F - t1F * tA[a]) == sp.zeros(16) for a in (1, 2, 3)))
kap = sp.Rational(2, 9)
for name, X in [('1 - k G^T', sp.eye(4) - kap * G.T), ('2 - G^T + 3 G^T^2', 2 * sp.eye(4) - G.T + 3 * G.T**2)]:
    rhoX = sp.kronecker_product(X, X.inv().T)
    comm = all((rhoX * words[(a, 0)] - words[(a, 0)] * rhoX) == sp.zeros(16) for a in (1, 2, 3)) and (rhoX * t1F - t1F * rhoX) == sp.zeros(16)
    v = X[0, :]; u = X.inv().T[3, :]
    sub = {thS: th1N - th2N, hS: hN}
    sub.update({vS[i]: v[i] for i in range(4)}); sub.update({uS[i]: u[i] for i in range(4)})
    cX = {k: sp.nsimplify(e.subs(sub)) for k, e in cc.items()}
    P = sum((cX[k] * words[k] for k in words), sp.zeros(16))
    print(name, ': rho(X) commutes with t:', comm, '; rho(X) == sum c_ak(X) t_a(th1-h) t1(th2)^k as operators:', (P - rhoX).applyfunc(sp.simplify) == sp.zeros(16),
          '; #nonzero c_ak:', sum(1 for x in cX.values() if x != 0))
