import sympy as sp, time
from chain import *
th1, th2, h = sp.symbols('theta1 theta2 h'); cs = sp.symbols('c1:5')
ch = Chain(['A', 'F'], [th1, th2], h, companion(cs))
t0 = time.time()
for a, u in [(2, th2), (3, th2), (4, th2), (4, th1 - h)]:
    M = ch.op_t(a, u).applyfunc(sp.expand)
    print('t_%d(%s) == 0 as operator:' % (a, u), M == sp.zeros(16), round(time.time() - t0, 1), flush=True)
# qKZ: t1(th1) = h rho1(G^T)(th12 - hK), t1(th2) = h rho2(G)(th12 - hK)  [bra action]; check
G = companion(cs)
K = sp.zeros(16)
for i in range(4):
    for j in range(4):
        K += sp.kronecker_product(E(i, j), E(i, j))
t1a = ch.op_t(1, th1); t1b = ch.op_t(1, th2)
lhsA = h * sp.kronecker_product(G.T, sp.eye(4)) * ((th1 - th2) * sp.eye(16) - h * K)
lhsB = h * sp.kronecker_product(sp.eye(4), G) * ((th1 - th2) * sp.eye(16) - h * K)
print('qKZ t1(th1) (bra-action matrix) == h rho1(G^T)(th12 - hK) [matrix order rho*R]:', (t1a - lhsA).applyfunc(sp.expand) == sp.zeros(16),
      ' or R*rho:', (t1a - h * ((th1 - th2) * sp.eye(16) - h * K) * sp.kronecker_product(G.T, sp.eye(4))).applyfunc(sp.expand) == sp.zeros(16))
print('qKZ t1(th2):', (t1b - lhsB).applyfunc(sp.expand) == sp.zeros(16), (t1b - h * ((th1 - th2) * sp.eye(16) - h * K) * sp.kronecker_product(sp.eye(4), G)).applyfunc(sp.expand) == sp.zeros(16))
