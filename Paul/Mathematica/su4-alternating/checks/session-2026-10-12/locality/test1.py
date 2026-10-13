import numpy as np, time
from lib import *
ch6 = build_chain(6)
G = ch6.G
print(np.abs(cof(G) - np.linalg.det(G)*np.linalg.inv(G).T).max())
rng = np.random.default_rng(1)
psi = rng.standard_normal(4096)
t=time.time(); out = ch6.apply(3, -6/5, psi); print('t3 time', time.time()-t)
t=time.time(); out = ch6.apply(1, -6/5, psi); print('t1 time', time.time()-t)
t=time.time(); out = ch6.apply(2, -6/5, psi); print('t2 time', time.time()-t)
us = {1:0.37, 2:-0.55, 3:-1.2}
def comm(a,b,ch,psi):
    x = ch.apply(b, 0.41, ch.apply(a, -0.37, psi))
    y = ch.apply(a, -0.37, ch.apply(b, 0.41, psi))
    return np.linalg.norm(x-y)/np.linalg.norm(x)
for ch,name in [(ch6,'6'),(build_chain(4),'4')]:
    p = rng.standard_normal(ch.dim)
    for a,b in [(1,2),(1,3),(2,3),(1,1),(2,2),(3,3)]:
        print(name, a,b, comm(a,b,ch,p))
