import sympy as sp, mpmath as mp, pickle, time
from L2sym import t, vac, h, d, Z, Gc
mp.mp.dps = 60
# symbolic operators in d (theta_A = d, theta_F = 0)
uu = sp.symbols('uu')
th = [d, 0]
t0 = time.time()
tsym = {a: t(a, uu, th) for a in range(1,5)}   # polynomial in uu, d
print('symbolic t_a(u) built', time.time()-t0, flush=True)
pickle.dump(tsym, open('L2tsym.pkl','wb'))
