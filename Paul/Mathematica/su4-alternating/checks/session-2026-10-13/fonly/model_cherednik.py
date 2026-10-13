# F-chain Cherednik model (python fallback used while the shared wolfbook kernel was blocked, 2026-10-13).
# Levels {0,1,2}^L (e4,e3,e2), N raises level (dropping level 3), row-vector action, Y_s = R_{s,s+1}..R_{s,L} N_s R_{s,1}..R_{s,s-1},
# R_{sf} = 1 - h/(th_s - th_f) P_{sf}.  Target tensor of f(N) = p - s N + N^2 on the vacuum; solve  sum_k c_k vac.Y^k = target.
import itertools, sympy as sp

def model(L, th, h):
    basis = list(itertools.product(range(3), repeat=L)); n = len(basis); idx = {t: i for i, t in enumerate(basis)}
    def P(i, j):
        M = sp.zeros(n, n)
        for q, t in enumerate(basis):
            u = list(t); u[i], u[j] = u[j], u[i]; M[q, idx[tuple(u)]] = 1
        return M
    def Nm(s):
        M = sp.zeros(n, n)
        for q, t in enumerate(basis):
            if t[s] < 2:
                u = list(t); u[s] += 1; M[q, idx[tuple(u)]] = 1
        return M
    R = lambda s, f: sp.eye(n) - (h / (th[s] - th[f])) * P(s, f)
    Y = []
    for s in range(L):
        M = sp.eye(n)
        for f in range(s + 1, L): M = M * R(s, f)
        M = M * Nm(s)
        for f in range(0, s): M = M * R(s, f)
        Y.append(M)
    return basis, idx, Y

def solve(L, th, h, p, s):
    basis, idx, Y = model(L, th, h)
    n = len(basis); vac = sp.zeros(1, n); vac[0, 0] = 1
    ks = list(itertools.product(range(3), repeat=L))
    rows = []
    for k in ks:
        v = vac
        for si, kk in enumerate(k):
            for _ in range(kk): v = v * Y[si]
        rows.append(v)
    A = sp.Matrix.vstack(*rows)                       # rows: words
    w = [p, -s, 1]
    tgt = sp.Matrix([[sp.prod([w[l] for l in t]) for t in basis]])
    c = (A.T).LUsolve(tgt.T)
    return ks, c

if __name__ == "__main__":
    h = sp.Rational(1, 2); p, s = sp.symbols('p s')
    th = [sp.Rational(1, 3), sp.Rational(-2, 5)]
    ks, c = solve(2, th, h, p, s)
    nu = th[0] - th[1]; g = nu**2 / (nu**2 - h**2)
    e1, e2 = sp.symbols('e1 e2'); x1, x2 = sp.symbols('x1 x2')
    F = sum(sp.simplify(ci) * x1**k[0] * x2**k[1] for ci, k in zip(c, ks))
    Phi = sp.expand((p - s*x1 + x1**2) * (p - s*x2 + x2**2))
    # Phi in e1,e2 then e2 -> g e2
    Pe = sp.expand(p**2 - s*p*e1 + p*(e1**2 - 2*e2) + s**2*e2 - s*e1*e2 + e2**2)
    Fpred = sp.expand(Pe.subs(e2, g*e2).subs({e1: x1 + x2, e2: x1*x2}))
    print('L=2: F - Phi(e1, g e2) =', sp.simplify(sp.expand(F - Fpred)))
