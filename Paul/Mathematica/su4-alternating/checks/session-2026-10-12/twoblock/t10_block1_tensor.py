import numpy as np
from core import *
from exact import *
R8 = res8(); p = list(R8.keys())[0]; Y = R8[p][1]
def idx(i1, i2, i3, i4): return ((i1 * 4 + i2) * 4 + i3) * 4 + i4
# T[i, j, a1, k2, k4] = Y[(i,j,0,3), (a1,k2,0,k4)]
T = np.zeros((4, 4, 4, 4, 4), dtype=np.int64)
for i in range(4):
    for j in range(4):
        for a1 in range(4):
            for k2 in range(4):
                for k4 in range(4):
                    T[i, j, a1, k2, k4] = Y[idx(i, j, 0, 3), idx(a1, k2, 0, k4)]
print('k4=3 slice zero:', not T[..., 3].any())
T = T[..., :3]
# flattenings
def fr(A, axes_row):
    axes_col = [a for a in range(A.ndim) if a not in axes_row]
    sh = A.shape
    Mx = np.transpose(A, axes_row + axes_col).reshape(int(np.prod([sh[a] for a in axes_row])), -1)
    return modrank(Mx, p)
print('rank (i,a1) | rest :', fr(T, [0, 2]))
print('rank (j,k2) | rest :', fr(T, [1, 3]))
print('rank k4 | rest     :', fr(T, [4]))
print('rank (i,a1,j,k2)|k4:', fr(T, [0, 2, 1, 3]))
print('rank (i,j) | rest  :', fr(T, [0, 1]))
print('rank a1 | rest     :', fr(T, [2]))
print('rank i | rest      :', fr(T, [0]))
# modulo the singlet: add delta_ij X[a1,k2,k4]; minimal rank? compute rank of (i,a1) flattening of the traceless projection
Tt = T.copy()
tr = sum(T[k, k] for k in range(4)) % p
inv4 = pow(4, -1, p)
for k in range(4):
    Tt[k, k] = (T[k, k] - tr * inv4) % p
print('traceless part: rank (i,a1)|rest', fr(Tt, [0, 2]), ' (j,k2)|rest', fr(Tt, [1, 3]), ' (i,j)|rest', fr(Tt, [0, 1]))
np.save('T_block1.npy', T)
