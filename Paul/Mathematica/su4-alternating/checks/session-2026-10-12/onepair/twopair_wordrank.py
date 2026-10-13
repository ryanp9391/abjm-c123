import numpy as np, pickle
from fastchain import *
d = pickle.load(open('twopair_W.pkl', 'rb')); labs = d['labs']; W4 = d['W']
rk = {x: rank_mod(W4[n].reshape(16, 16)) for n, x in enumerate(labs)}
print('rank of each 4-site word as pair1 x pair2 matrix; rows (a1,k2), cols (a3,k4):')
for a1 in range(4):
    for k2 in range(4):
        print((a1, k2), [rk[(a1, k2, a3, k4)] for a3 in range(4) for k4 in range(4)])
# is <0|t_a1(th1-h) t1(th2)^k2 = <(a1,k2)|_12 (x) <0_34| (the local 2-site word)?
T = pickle.load(open('twopair_T.pkl', 'rb'))
W12, W34 = T['W12'], T['W34']
vac34 = np.zeros(16, dtype=object); vac34[3] = 1
vac12 = np.zeros(16, dtype=object); vac12[3] = 1
ok1 = [x[:2] for n, x in enumerate(labs) if x[2:] == (0, 0) and np.array_equal(W4[n] % P, np.kron(W12[4 * x[0] + x[1]], vac34) % P)]
ok2 = [x[2:] for n, x in enumerate(labs) if x[:2] == (0, 0) and np.array_equal(W4[n] % P, np.kron(vac12, W34[4 * x[2] + x[3]]) % P)]
print('pair-1 words equal to local 2-site word (x) <0_34| :', ok1)
print('pair-2 words equal to <0_12| (x) local 2-site word :', ok2)
