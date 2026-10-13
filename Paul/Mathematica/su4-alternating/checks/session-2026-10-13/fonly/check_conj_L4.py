# Check at L = 3 (other theta) and L = 4: F = B_{z2} B_{z3}, B_z = sum_S (-z)^{L-|S|} g_S Y^S, g_S = prod_{s<t in S} nu^2/(nu^2-h^2)
import itertools, sympy as sp, sys
from sympy import Rational as Q
def run(L, th, h, z2, z3):
    p, s = z2*z3, z2 + z3
    basis = list(itertools.product(range(3), repeat=L)); n = len(basis); idx = {t: i for i, t in enumerate(basis)}
    def P(i, j):
        M = sp.zeros(n, n)
        for q, t in enumerate(basis):
            u = list(t); u[i], u[j] = u[j], u[i]; M[q, idx[tuple(u)]] = 1
        return M
    def Nm(sx):
        M = sp.zeros(n, n)
        for q, t in enumerate(basis):
            if t[sx] < 2:
                u = list(t); u[sx] += 1; M[q, idx[tuple(u)]] = 1
        return M
    R = lambda i, j: sp.eye(n) - (h/(th[i] - th[j]))*P(i, j)
    Y = []
    for i in range(L):
        M = sp.eye(n)
        for f in range(i+1, L): M = M*R(i, f)
        M = M*Nm(i)
        for f in range(0, i): M = M*R(i, f)
        Y.append(M)
    vac = sp.zeros(1, n); vac[0, 0] = 1
    def g(i, j): nu = th[i]-th[j]; return nu**2/(nu**2 - h**2)
    def Bop(z):
        tot = sp.zeros(n, n)
        for r in range(L+1):
            for S in itertools.combinations(range(L), r):
                gS = sp.prod([g(i, j) for i, j in itertools.combinations(S, 2)]) if len(S) > 1 else 1
                T = sp.eye(n)
                for i in S: T = T*Y[i]
                tot += (-z)**(L-r)*gS*T
        return tot
    lhs = vac*Bop(z2)*Bop(z3)
    w = [p, -s, 1]
    tgt = sp.Matrix([[sp.prod([w[l] for l in t]) for t in basis]])
    return (lhs - tgt).norm() == 0
print('L=3 other theta:', run(3, [Q(2,7), Q(-3,5), Q(8,11)], Q(1,3), Q(5,4), Q(2,3)), flush=True)
print('L=4:', run(4, [Q(1,3), Q(-2,5), Q(3,7), Q(-5,11)], Q(1,2), Q(3,2), Q(-2,3)), flush=True)
