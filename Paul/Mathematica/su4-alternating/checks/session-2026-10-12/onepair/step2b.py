import sympy as sp, pickle
U = sp.Symbol('U')
hN = sp.Rational(1, 2); th1N, th2N = sp.Rational(2, 7), sp.Rational(-3, 5)
zN = [sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(4, 7)]
d = pickle.load(open('step2.pkl', 'rb')); res = d['res']
V4 = lambda x: (x - th1N) * (x - th2N - hN)
def baxter_q(taus, zj, D=4):
    cs = sp.symbols('q0:%d' % (D + 1))
    q = lambda x: sum(cs[n] * x**n for n in range(D + 1))
    tt = {0: sp.Integer(1)}; tt.update({a: taus[a].as_expr() for a in range(1, 5)})
    expr = sum((-1)**a * tt[a] * zj**(4 - a) * q(U + hN - a * hN) * sp.prod([V4(U - b * hN) for b in range(a, 4)]) for a in range(5))
    eqs = sp.Poly(sp.expand(expr), U).all_coeffs()
    sol = sp.linsolve(eqs, cs)
    return sol, cs
for key in list(res)[:12]:
    out = []
    for j in range(4):
        sol, cs = baxter_q(res[key], zN[j])
        s = list(sol)[0]
        free = set().union(*[sp.sympify(x).free_symbols for x in s])
        poly = sp.expand(sum(s[n] * U**n for n in range(len(cs))))
        out.append((j + 1, sp.Poly(poly, U).degree() if poly != 0 else None, len(free)))
    print(key, out)
