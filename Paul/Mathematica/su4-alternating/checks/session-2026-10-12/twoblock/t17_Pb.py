import numpy as np, pickle
from core import *
pieces = pickle.load(open('phys_pieces_modp.pkl', 'rb')); p = P1
for k, c in pieces.items():
    c4 = c.reshape(4, 4, 4, 4)[:, :, 0, :3]  # a1, k2, k4
    print(k, [ (a1, modrank(c4[a1], p), [(k2, k4) for k2 in range(4) for k4 in range(3) if c4[a1, k2, k4]]) for a1 in range(4)])
    # site-2 x (site1,site4) and site-4 x (site1, site2) flattening ranks
    print('   rank k2|(a1,k4):', modrank(np.transpose(c4, (1, 0, 2)).reshape(4, 12), p), ' rank k4|(a1,k2):', modrank(np.transpose(c4, (2, 0, 1)).reshape(3, 16), p), ' rank a1|(k2,k4):', modrank(c4.reshape(4, 12), p))
