import numpy as np, pickle
from collections import Counter
from core import *
res = pickle.load(open('dual_mod.pkl', 'rb')); M, Y = res[P1]
nz = np.argwhere(Y != 0)
deg = lambda idx: (idx[0]) + (3 - idx[1]) + (idx[2]) + (3 - idx[3])   # 0-based: A: i-1, F: 4-j
viol = [(wordsIdx[r], wordsIdx[c]) for r, c in nz if deg(wordsIdx[r]) != sum(wordsIdx[c])]
print('Y entries:', len(nz), ' grading violations:', len(viol), viol[:5])
nzM = np.argwhere(M != 0)
violM = [(wordsIdx[r], wordsIdx[c]) for r, c in nzM if sum(wordsIdx[r]) != deg(wordsIdx[c])]
print('M entries:', len(nzM), ' grading violations:', len(violM))
# block sizes by degree
cnt_lab = Counter(sum(x) for x in wordsIdx); cnt_comp = Counter(deg(x) for x in wordsIdx)
print('labels per degree :', sorted(cnt_lab.items()))
print('components per deg:', sorted(cnt_comp.items()))
# is M block-triangular within degree w.r.t. some order? nnz per degree block
for d in range(13):
    rows = [i for i, x in enumerate(wordsIdx) if sum(x) == d]; cols = [i for i, x in enumerate(wordsIdx) if deg(x) == d]
    blk = M[np.ix_(rows, cols)]
    print(' degree', d, 'block', blk.shape, 'nnz', np.count_nonzero(blk), 'nnz inverse', np.count_nonzero(Y[np.ix_(cols, rows)]))
