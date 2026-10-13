import numpy as np, pickle
from fastchain import *
d = pickle.load(open('twopair_W.pkl', 'rb')); labs = d['labs']; W = d['W']
Winv = pickle.load(open('twopair_Winv.pkl', 'rb'))
from collections import Counter
rk = {}
for n, x in enumerate(labs):
    Mx = Winv[:, n].reshape(16, 16)   # rows (i1 j1), cols (i2 j2)
    rk[x] = rank_mod(Mx)
print('distribution of ranks of the pair1 x pair2 bilinear form of each coefficient c_x:', Counter(rk.values()))
print('words with rank 1:', [x for x in labs if rk[x] == 1])
print('ranks table rows (a1,k2), cols (a3,k4):')
for a1 in range(4):
    for k2 in range(4):
        print((a1, k2), [rk[(a1, k2, a3, k4)] for a3 in range(4) for k4 in range(4)])
