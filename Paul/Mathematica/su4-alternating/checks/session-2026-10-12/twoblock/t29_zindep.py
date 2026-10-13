import numpy as np
from core import *
R = Ring('mod', P1)
M1, _ = sov_words(Chain(R))
M2, _ = sov_words(Chain(R, z=[Fr(7, 3), Fr(-2, 5), Fr(5, 6), Fr(9, 7)]))
print('words z-independent (exact mod p):', np.array_equal(M1, M2))
