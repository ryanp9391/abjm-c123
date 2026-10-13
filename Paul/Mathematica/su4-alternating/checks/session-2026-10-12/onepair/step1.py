import sympy as sp, pickle, time
from chain import *
th1, th2, h = sp.symbols('theta1 theta2 h')
cs = sp.symbols('c1:5')
G = companion(cs)
ch = Chain(['A', 'F'], [th1, th2], h, G)
vac = sp.zeros(1, 16); vac[0, 0*4 + 3] = 1   # e1 (x) e4
t0 = time.time()
t1th2 = ch.op_t(1, th2)
print('t1(th2) built', time.time() - t0)
rows = {}
for a in range(4):
    br = ch.bra_t(vac, a, th1 - h)
    for k in range(4):
        rows[(a, k)] = br.applyfunc(sp.factor)
        br = (br * t1th2).applyfunc(sp.expand)
print('basis built', time.time() - t0)
W = sp.Matrix.vstack(*[rows[(a, k)] for a in range(4) for k in range(4)])
print('free symbols of W:', W.free_symbols)
pickle.dump({'W': W, 't1th2': t1th2}, open('step1.pkl', 'wb'))
