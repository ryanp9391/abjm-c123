import pickle, itertools, sympy as sp
d = pickle.load(open('L3_coeffs_th.pkl', 'rb'))
p, s, z2, z3 = sp.symbols('p s z2 z3')
h = sp.Rational(1, 2); th = [sp.Rational(1, 3), sp.Rational(-2, 5), sp.Rational(3, 7)]
L = 3; x = sp.symbols('x1 x2 x3')
def g(i, j):
    nu = th[i] - th[j]; return nu**2/(nu**2 - h**2)
def B(z):
    tot = 0
    for r in range(L + 1):
        for S in itertools.combinations(range(L), r):
            gS = sp.prod([g(i, j) for i, j in itertools.combinations(S, 2)]) if len(S) > 1 else 1
            tot += (-z)**(L - r) * gS * sp.prod([x[i] for i in S])
    return tot
F = sp.expand(B(z2) * B(z3))
F = sp.expand(F.subs({z2: (s + sp.sqrt(s**2 - 4*p))/2, z3: (s - sp.sqrt(s**2 - 4*p))/2}))
F = sp.simplify(F)
bad = 0
for k in itertools.product(range(3), repeat=3):
    ck = sp.Poly(F, *x).coeff_monomial(x[0]**k[0] * x[1]**k[1] * x[2]**k[2])
    ex = sp.sympify(d[str(k)])
    if sp.simplify(ck - ex) != 0: bad += 1; print('mismatch', k, sp.simplify(ck), ex)
print('mismatches:', bad, 'of 27')
