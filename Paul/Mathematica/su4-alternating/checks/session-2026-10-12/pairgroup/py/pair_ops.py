# 1-pair chain (A,F) at the ABJM point theta_A = theta_F + 2h: operators X1 = t1(thF), X2 = t2(thF+h), X3 = t3hat(thF+h)
# question: is rho_pair(1 - kappa X) in span{1, X1, X2, X3, ...} ?
import numpy as np
from pgbase import *
rng = np.random.default_rng(3)
h = 0.5; thF = -0.6; thA = thF + 2*h
def opmat(ch, k, u, drop=None):
    """full matrix of t_k(u) on the chain (bras: psi.M), optionally divided by vanishing factor"""
    N = ch.N; dim = 4**N; M = np.zeros((dim, dim))
    for i in range(dim):
        e = np.zeros(dim); e[i] = 1
        M[i] = ch.apply_t(k, u, e.reshape((4,)*N)).reshape(-1)
    return M
def t3hat_pair(ch, u):
    # t3 without the prefactor (u-thF)(u-thF-h) of the F site of pair 1, keep A prefactor
    return None
G = rng.normal(size=(4,4)); Gi = np.linalg.inv(G)
ch = Chain(['A','F'], [thA, thF], h, G.astype(object))
X1 = opmat(ch, 1, thF)
X2 = opmat(ch, 2, thF + h)
# t3hat: evaluate t3 at thF+h+eps / eps  (numerical limit) -- better: exact formula without g_F. Do via finite eps symmetric
eps = 1e-6
X3 = (opmat(ch, 3, thF + h + eps) - opmat(ch, 3, thF + h - eps))/(2*eps)
X4 = opmat(ch, 4, 0.123)
print('X4 scalar?', np.abs(X4 - X4[0,0]*np.eye(16)).max())
def rhopair(g): return np.kron(g, np.linalg.inv(g).T)
basis = [np.eye(16), X1, X2, X3]
names = ['1','X1','X2','X3']
B = np.array([b.reshape(-1) for b in basis]).T
print('rank of {1,X1,X2,X3}:', np.linalg.matrix_rank(B, 1e-8))
for Xname, X in (('G',G),('G^T',G.T),('G^-1',Gi),('G^-T',Gi.T)):
    for kap in (0.3, -0.7):
        tgt = rhopair(np.eye(4) - kap*X).reshape(-1)
        c, *_ = np.linalg.lstsq(B, tgt, rcond=None)
        print(Xname, kap, 'resid', np.linalg.norm(B@c - tgt)/np.linalg.norm(tgt))
# products
prods = {'X1X1':X1@X1, 'X1X2':X1@X2, 'X1X3':X1@X3, 'X2X3':X2@X3, 'X2X2':X2@X2, 'X3X3':X3@X3}

print('--- on BPS states ---')
def bps():
    v = rng.normal(size=4); u = rng.normal(size=4); u = u - (v@u)/(v@v)*v; return v, u
v, u = bps()
psi = np.kron(v, u)
for setname, ops in (('1,X1,X2,X3', [np.eye(16), X1, X2, X3]),
                     ('+X1X3', [np.eye(16), X1, X2, X3, X1@X3]),
                     ('+X1X3,X1X1,X1X2,X2X3', [np.eye(16), X1, X2, X3, X1@X3, X1@X1, X1@X2, X2@X3])):
    B = np.array([psi@o for o in ops]).T
    print(setname, 'rank', np.linalg.matrix_rank(B, 1e-9))
    for Xname, X in (('G',G),('G^T',G.T),('G^-1',Gi),('G^-T',Gi.T)):
        for kap in (0.3, -0.7):
            tgt = psi @ rhopair(np.eye(4) - kap*X)
            c, *_ = np.linalg.lstsq(B, tgt, rcond=None)
            print('  ', Xname, kap, 'resid %.2e'%(np.linalg.norm(B@c - tgt)/np.linalg.norm(tgt)), np.round(c,5))
