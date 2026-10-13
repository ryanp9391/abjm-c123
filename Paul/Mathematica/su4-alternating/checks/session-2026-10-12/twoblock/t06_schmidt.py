import numpy as np, pickle
from core import *
res = pickle.load(open('dual_mod.pkl', 'rb'))
p = P1; M, Y = res[p]
Y8 = Y.reshape(4,4,4,4, 4,4,4,4)   # i1 i2 i3 i4 | a1 k2 a3 k4
M8 = M.reshape(4,4,4,4, 4,4,4,4)   # a1 k2 a3 k4 | i1 i2 i3 i4
def schmidt(T8, groupA, name):
    # T8 axes: 0..3 product-site indices, 4..7 label slots (slot s <-> site s)
    axesA = [s for s in groupA] + [4 + s for s in groupA]
    axesB = [s for s in range(4) if s not in groupA] + [4 + s for s in range(4) if s not in groupA]
    Mx = np.transpose(T8, axesA + axesB).reshape(4 ** (2 * len(groupA)), -1)
    print(name, 'Schmidt rank cut', [s + 1 for s in groupA], '|', [s + 1 for s in range(4) if s not in groupA], ':', modrank(Mx, p))
for g in ([0, 1], [0], [3], [0, 1, 2], [1], [2], [0, 2]):
    schmidt(Y8, g, 'Y=M^-1')
M8t = np.transpose(M8, (4,5,6,7,0,1,2,3))
for g in ([0, 1], [0], [3], [0,1,2]):
    schmidt(M8t, g, 'M      ')
