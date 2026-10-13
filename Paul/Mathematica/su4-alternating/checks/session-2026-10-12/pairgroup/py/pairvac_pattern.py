# Task 3: <e1 e4| (L_A L_F)_{jk}(u) = (u-thA)(u-thF) d_jk <e1e4| - h(u-thA) d_k4 <e1 e_j| + h(u-thF) d_j1 <e_k e4| - h^2 d_j1 d_k4 <S|
import sympy as sp, numpy as np
from pgbase import lax_k
u, h, tA, tF = sp.symbols('u h thA thF')
LA = lax_k('A', 1, u, tA, h); LF = lax_k('F', 1, u, tF, h)
vac = sp.zeros(1, 16); vac[0, 0*4+3] = 1
def e(i): v = sp.zeros(1, 4); v[0, i] = 1; return v
def kron(a, b): return sp.Matrix(1, 16, lambda _, n: a[0, n//4]*b[0, n % 4])
S = sum((kron(e(i), e(i)) for i in range(4)), sp.zeros(1, 16))
ok = True; zeros = {}
for j in range(4):
    for k in range(4):
        M = sp.zeros(16, 16)
        for m in range(4):
            A = sp.Matrix(LA[j, m].tolist()); F = sp.Matrix(LF[m, k].tolist())
            M += sp.kronecker_product(A, F)
        lhs = (vac*M).applyfunc(sp.expand)
        rhs = ((u-tA)*(u-tF)*(1 if j == k else 0)*vac - h*(u-tA)*(1 if k == 3 else 0)*kron(e(0), e(j))
               + h*(u-tF)*(1 if j == 0 else 0)*kron(e(k), e(3)) - h**2*(1 if (j == 0 and k == 3) else 0)*S).applyfunc(sp.expand)
        ok &= (lhs - rhs).is_zero_matrix
        zeros[(j+1, k+1)] = ('0' if lhs.is_zero_matrix else ('0 at u=thA' if lhs.subs(u, tA).is_zero_matrix else ('0 at u=thF' if lhs.subs(u, tF).is_zero_matrix else 'nonzero')))
print('formula holds for all (j,k):', ok)
for jk, s in zeros.items(): print(jk, s)
