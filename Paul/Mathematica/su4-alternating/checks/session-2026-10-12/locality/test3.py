import numpy as np
from lib import *
def site_apply(psi, s, M, n):
    return np.einsum('xay,ab->xby', psi.reshape(4**s,4,4**(n-s-1)), M).reshape(-1)
n=4
ch = build_chain(n)
th = [float(x) for x in ch.th]; h=0.5
c1inv = h
for s in range(1,n):
    c1inv *= (th[0]-th[s]) if ch.types[s]=='F' else (th[0]-th[s]+h)
c1=1/c1inv
Li = np.diag([1/float(x) for x in Z4])
rng=np.random.default_rng(5)
def test(psi):
    lhs = c1*site_apply(psi,0,Li,n); rhs = ch.apply(1, th[0], lhs); return np.linalg.norm(rhs-psi)/np.linalg.norm(psi)
for a in (1/3,-2/3,0.7):
    print('W_a^L', a, test(W_state(n,a)))
# product with all sites supported on {1,4} random
vs=[]
for s in range(n):
    v=np.zeros(4); v[0]=rng.standard_normal(); v[3]=rng.standard_normal(); vs.append(v)
print('random 1,4 support', test(kron_all(vs)))
# BPS supported {1,4}: v=(1,0,0,a), u=(-b,0,0,1) v.u=-b+a
vs=[]
for p in range(2):
    a=rng.standard_normal(); v=np.array([1,0,0,a]); u=np.array([-a*rng.standard_normal()+0,0,0,1.0]); 
    vs += [v,u]
print('(1,a)x(c,1) support, not BPS', test(kron_all(vs)))
# pair 1 W_a', other pairs W_a
print('W_{a1} x W_a', test(W_state(n,1/3,a1=-2/3)))
# only A-site1 random, rest W
v,u = site_vec(1/3)
for typ in ['site1 random']:
    vs=[rng.standard_normal(4), u] + [v,u]
    print(typ, test(kron_all(vs)))
# BPS with  v general, pair 2 W
vv = rng.standard_normal(4); uu = rng.standard_normal(4); uu -= vv*(vv@uu)/(vv@vv)
print('BPS pair 1 random, pair 2 W_a', test(kron_all([vv,uu,v,u])))
print('BPS pair 2 random, pair 1 W_a', test(kron_all([v,u,vv,uu])))
