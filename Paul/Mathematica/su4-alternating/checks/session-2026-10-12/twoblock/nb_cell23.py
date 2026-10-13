import numpy as np, pickle, math
from core import *
d = np.load('nb_states_std.npz'); V = d['V']; vac = d['vac']; A1 = d['A1']; r2 = d['r2']; r4 = d['r4']
res = pickle.load(open('phys_pieces_exact.pkl', 'rb'))
b = Fr(3, 7)
c = {}
for k, vals in res.items():
    for x, v in vals.items(): c[x] = c.get(x, 0) + v * b ** k
z = z0
lam = [Fr(1) / math.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
Vf = np.array([[float(x) for x in row] for row in Vand(z)]); Vi = np.linalg.inv(Vf)
Bc = np.kron(np.kron(np.eye(4)[0] + float(b * lam[0]) * Vf[:, 3], np.eye(4)[3] - float(b) * Vi[0]), np.kron(np.eye(4)[0], np.eye(4)[3]))
Kdir = (Bc @ V) / (vac @ V)
Ksov = sum(float(v) * A1[x[0]] * r2 ** x[1] * r4 ** x[3] for x, v in c.items())
print('nonzero c:', len(c), ' a3 all zero:', all(x[2] == 0 for x in c))
print('max rel diff K_sov vs K_dir over 256 states (float64): %.2e' % (np.abs(Ksov - Kdir) / np.abs(Kdir)).max())
