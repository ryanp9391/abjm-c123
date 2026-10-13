import itertools, pickle
import numpy as np
exec(open('ketC.py').read().split("for rowset in")[0])
for rowset in ([0, 1, 2], [1, 2, 3]):
    for colsgn in (-1, 1):
        for argsgn in (-1, 1):
            C1 = Cop(0.37, colsgn, argsgn, rowset)
            print(rowset, colsgn, argsgn, '|C| =', np.abs(C1).max())
