"""Physical-frame images (exact sympy): (i) the would-be common rotation 1 - kappa Gc^T on <0| pair 1;
(ii) 2603-type dressed references (A: (zeta - Gc^T)^{-1}, F: (zeta' - Gc)^{-1}) as zeta -> 1 with z = 1 + eps phi."""
import sympy as sp
e, k, c, xi, xi2 = sp.symbols('epsilon kappa c xi xi2')
phi = [sp.Rational(9, 10), sp.Rational(7, 20), sp.Rational(-3, 10), sp.Rational(-11, 10)]
z = [1 + e * p for p in phi]
lam = [1 / sp.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
r = [sp.simplify(l / lam[0]) for l in lam]
print('r_j (eps-independent):', r)
n = sp.Matrix([[1, 1, 1, r[3]], [0, 1, 0, -r[1]], [0, 0, 1, -r[2]], [0, 0, 0, 1]]); Dp = sp.diag(1, r[1], r[2], 1)
g0 = n * Dp
V = sp.Matrix(4, 4, lambda i, j: z[j] ** (3 - i))
Lam = sp.diag(*z)
Gc = sp.simplify(V * Lam * V.inv())
physA = lambda vc: sp.simplify(vc * V.inv().T * g0.inv() / lam[0])
physF = lambda uc: sp.simplify(uc * V * g0.T)
e1 = sp.Matrix([[1, 0, 0, 0]]); e4 = sp.Matrix([[0, 0, 0, 1]])
# (i) common rotation in companion frame D_c = 1 - kappa Gc^T  <->  diag 1 - kappa Lam  <->  physical n'(1-kappa Lam)n'^-1
M = sp.simplify(n * (sp.eye(4) - k * Lam) * n.inv())
print('(n(1-kLam)n^-1)_14 =', sp.simplify(M[0, 3]))
vA = physA(e1 * (sp.eye(4) - k * Gc.T)); uF = physF(e4 * (sp.eye(4) - k * Gc.T).inv().T)
print('common rotation, physical A vector:', [sp.factor(x) for x in vA])
print('common rotation, physical F vector (u3)_1/(u3)_4 =', sp.factor(uF[0] / uF[3]))
# (ii) dressed references, zeta = 1 + eps xi
zeta = 1 + e * xi
vA = physA(e1 * (zeta * sp.eye(4) - Gc.T).inv())
lead = [sp.limit(x * e ** 1 / vA[3] * vA[3] / e, e, 0) if False else sp.series(x, e, 0, 1).removeO() for x in vA]
print('A dressed (zeta - Gc^T)^-1, physical v (normalised by v_4), eps -> 0:', [sp.factor(sp.limit(x / vA[3], e, 0)) for x in vA])
zeta2 = 1 + e * xi2
uF = physF(e4 * (zeta2 * sp.eye(4) - Gc).inv())
print('F dressed (zeta\' - Gc)^-1, physical u (normalised by u_4), eps -> 0:', [sp.factor(sp.limit(x / uF[3], e, 0)) for x in uF])
