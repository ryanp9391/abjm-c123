import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; R = X.R
I = R.arr(np.eye(4, dtype=int))
G = X.Gc; GT = X.GcT
e1, e4 = np.array(e(1)), np.array(e(4))
def lin(k, A): return (I - (X.r(k) * A) % p) % p
def show(name, vs):
    c = X.coeffs(X.bra(vs)); c4 = c.reshape(4, 4, 4, 4)
    sup = [wordsIdx[i] for i in np.nonzero(c)[0]]
    print(f'{name:55s}', X.ranks(c))
    print('     support:', sup if len(sup) <= 40 else str(len(sup)) + ' labels')
k1, k2, k3, k4 = Fr(3, 7), Fr(-5, 11), Fr(2, 9), Fr(7, 5)
show('site1: e1(1-k G^T)', [X.vm(e1, lin(k1, GT)), e4, e1, e4])
show('site2: e4(1-k G)', [e1, X.vm(e4, lin(k2, G)), e1, e4])
show('site2: e4(1-k G)^-1', [e1, X.vm(e4, modinv(lin(k2, G), p)), e1, e4])
show('site3: e1(1-k G^T)', [e1, e4, X.vm(e1, lin(k3, GT)), e4])
show('site4: e4(1-k G)', [e1, e4, e1, X.vm(e4, lin(k4, G))])
show('sites1,2: e1(1-k1 G^T) x e4(1-k2 G)', [X.vm(e1, lin(k1, GT)), X.vm(e4, lin(k2, G)), e1, e4])
show('sites3,4: e1(1-k3 G^T) x e4(1-k4 G)', [e1, e4, X.vm(e1, lin(k3, GT)), X.vm(e4, lin(k4, G))])
show('pair1 common g=1-kG^T: v g x u g^-T', [X.vm(e1, lin(k1, GT)), X.vm(e4, modinv(lin(k1, G), p)), e1, e4])
show('pair2 common g: ', [e1, e4, X.vm(e1, lin(k3, GT)), X.vm(e4, modinv(lin(k3, G), p))])
show('pair1 common g1, pair2 common g2', [X.vm(e1, lin(k1, GT)), X.vm(e4, modinv(lin(k1, G), p)), X.vm(e1, lin(k3, GT)), X.vm(e4, modinv(lin(k3, G), p))])
show('pair1 common g, pair2 common g (global)', [X.vm(e1, lin(k1, GT)), X.vm(e4, modinv(lin(k1, G), p)), X.vm(e1, lin(k1, GT)), X.vm(e4, modinv(lin(k1, G), p))])
show('all four sites linear', [X.vm(e1, lin(k1, GT)), X.vm(e4, lin(k2, G)), X.vm(e1, lin(k3, GT)), X.vm(e4, lin(k4, G))])
