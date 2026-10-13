import numpy as np, pickle, random
from fastchain import *
d = pickle.load(open('twopair_W.pkl', 'rb')); labs = d['labs']
Winv = pickle.load(open('twopair_Winv.pkl', 'rb'))
def coeffs(v1, u1, v2, u2):
    x = np.kron(np.kron(np.kron(np.array(v1, dtype=object), np.array(u1, dtype=object)), np.array(v2, dtype=object)), np.array(u2, dtype=object)) % P
    c = x.dot(Winv) % P
    C = np.zeros((16, 16), dtype=object)
    for n, (a1, k2, a3, k4) in enumerate(labs): C[4 * a1 + k2, 4 * a3 + k4] = c[n]
    return C
random.seed(5)
rnd = lambda: [random.randrange(P) for _ in range(4)]
e = lambda i: [1 if k == i else 0 for k in range(4)]
def show(C, title):
    print(title)
    print('      cols (a3,k4) = ' + ' '.join('%d%d' % (a, k) for a in range(4) for k in range(4)))
    for r in range(16):
        print('  (%d,%d) ' % divmod(r, 4), ' '.join(' X' if C[r, c] % P else ' .' for c in range(16)))
show(coeffs(rnd(), rnd(), e(0), e(3)), 'pair 1 generic, pair 2 vacuum')
show(coeffs(e(0), e(3), rnd(), rnd()), 'pair 1 vacuum, pair 2 generic')
show(coeffs(e(0), e(3), rnd(), e(3)), 'pair 1 vacuum, pair 2 = (v2, e4) [A-site 3 rotated only]')
show(coeffs(e(0), e(3), e(0), rnd()), 'pair 1 vacuum, pair 2 = (e1, u2) [F-site 4 rotated only]')
