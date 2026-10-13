import sympy as sp, pickle
U = sp.Symbol('U')
zN = [sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(4, 7)]
d = pickle.load(open('step2.pkl', 'rb')); tops = d['tops']
Vv = lambda z: sp.Matrix([z**3, z**2, z, 1])
def pvec(a):
    others = [zN[m] for m in range(4) if m != a]
    return sp.Matrix(sp.Poly(sp.prod([U - zz for zz in others]), U).all_coeffs())
B = sp.Matrix.hstack(*[sp.kronecker_product(pvec(i), Vv(zN[i])) for i in range(4)])
u0 = sp.Rational(3, 11)
M1 = tops[1].subs(U, u0)
# restrict: M1*B = B*R
R = (B.T * B).inv() * B.T * (M1 * B)
print('closure:', (M1 * B - B * R).applyfunc(sp.simplify) == sp.zeros(16, 4))
lam = sp.Symbol('lam')
cp = sp.factor(R.charpoly(lam).as_expr())
print('char poly of t1(3/11) on zero-weight space:', cp)
# symbolic u: is R(u) = A + u*B..? check commuting family
Ru = (B.T * B).inv() * B.T * (tops[1] * B)
Ru = Ru.applyfunc(sp.factor)
print(Ru)
