# all 16 eigenstates of the A F chain at the rational point: exact (12 product) and mpmath 50 digits (4 mixed)
import sympy as sp, pickle, mpmath as mp
mp.mp.dps = 50
U = sp.Symbol('U')
hN = sp.Rational(1, 2); th1N, th2N = sp.Rational(2, 7), sp.Rational(-3, 5)
zN = [sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(4, 7)]
d = pickle.load(open('step2.pkl', 'rb')); tops = d['tops']
Vv = lambda z: sp.Matrix([z**3, z**2, z, 1])
def pvec(a):
    others = [zN[m] for m in range(4) if m != a]
    return sp.Matrix(sp.Poly(sp.prod([U - zz for zz in others]), U).all_coeffs())
def tau_polys_exact(psi):
    out = {}
    n = next(i for i in range(16) if psi[i] != 0)
    for a in range(1, 5):
        w = tops[a] * psi
        out[a] = [sp.Rational(c) for c in sp.Poly(sp.expand(w[n] / psi[n]), U).all_coeffs()[::-1]]  # ascending
    return out
states = []
for a in range(4):
    for b in range(4):
        if a != b:
            psi = sp.kronecker_product(pvec(a), Vv(zN[b]))
            states.append({'label': 'P%d%d' % (a + 1, b + 1), 'exact': True, 'ket': psi, 'tau': tau_polys_exact(psi)})
# mixed: basis B_i = p_i (x) V_i ; t1 restricted = f(U) + C
B = sp.Matrix.hstack(*[sp.kronecker_product(pvec(i), Vv(zN[i])) for i in range(4)])
Ru = ((B.T * B).inv() * B.T * (tops[1] * B)).applyfunc(sp.expand)
C = Ru.applyfunc(lambda e: e.subs(U, 0)) - Ru[0, 0].subs(U, 0) * sp.eye(4)
Cm = mp.matrix([[mp.mpf(sp.Rational(C[i, j]).p) / C[i, j].q if C[i, j] != 0 else 0 for j in range(4)] for i in range(4)])
E, ER = mp.eig(Cm)
Bm = mp.matrix([[mp.mpf(B[i, j].p) / B[i, j].q for j in range(4)] for i in range(16)])
for n in range(4):
    x = ER[:, n]
    psi = Bm * x
    # tau polys numerically: evaluate t_a(U) on psi at 9 points and fit (degree <= 8)
    tau = {}
    for a in range(1, 5):
        Ta = tops[a]
        deg = sp.Poly(Ta[3, 3] + Ta[0, 0] + Ta[5, 5], U).degree()
        pts = [mp.mpf(k) / 7 for k in range(-4, 5)]
        vals = []
        for p in pts:
            Mp = Ta.subs(U, sp.Rational(int(p * 7), 7))
            Mm = mp.matrix([[mp.mpf(Mp[i, j].p) / Mp[i, j].q for j in range(16)] for i in range(16)])
            w = Mm * psi
            k0 = max(range(16), key=lambda i: abs(psi[i]))
            vals.append(w[k0] / psi[k0])
        # polynomial fit degree 8 (exact interpolation)
        A = mp.matrix([[p**m for m in range(9)] for p in pts])
        coeffs = mp.lu_solve(A, mp.matrix(vals))
        tau[a] = [coeffs[m] for m in range(9)]
        # check eigen at one more point
    states.append({'label': 'M%d' % (n + 1), 'exact': False, 'ket': [psi[i] for i in range(16)], 'tau': tau, 'C_eig': E[n]})
pickle.dump(states, open('states.pkl', 'wb'))
for s in states:
    print(s['label'], 'tau1 coeffs:', [mp.nstr(mp.mpf(c) if not isinstance(c, sp.Basic) else mp.mpf(c.p) / c.q, 8) for c in s['tau'][1][:3]])
