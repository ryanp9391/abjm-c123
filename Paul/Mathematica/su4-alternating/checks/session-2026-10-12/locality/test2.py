import numpy as np
from lib import *
def site_apply(psi, s, M, n):
    return np.einsum('xay,ab->xby', psi.reshape(4**s,4,4**(n-s-1)), M).reshape(-1)
for n in (4,6):
    ch = build_chain(n)
    th = [float(x) for x in ch.th]; h=0.5
    c1inv = h
    for s in range(1,n):
        if ch.types[s]=='F': c1inv *= th[0]-th[s]
        else: c1inv *= (th[0]-th[s]+h)
    c1 = 1/c1inv
    print(n, 'c1', c1, 1/c1 if n==4 else '')
    Li = np.diag([1/float(x) for x in Z4])
    rng = np.random.default_rng(3)
    cases = {'W^L': W_state(n, 1/3)}
    vs = [rng.standard_normal(4) for _ in range(n)]
    cases['random product'] = kron_all(vs)
    # BPS product: pairs with v.u=0 random
    vs=[]
    for p in range(n//2):
        v = rng.standard_normal(4); u = rng.standard_normal(4); u -= v*(v@u)/(v@v)
        vs += [v,u]
    cases['random BPS product'] = kron_all(vs)
    for name, psi in cases.items():
        lhs = c1*site_apply(psi,0,Li,n)
        rhs = ch.apply(1, th[0], lhs)
        print(n, name, np.linalg.norm(rhs-psi)/np.linalg.norm(psi))
