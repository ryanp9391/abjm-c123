import sys
from pipeline import *
L = int(sys.argv[1]); tag = sys.argv[2]; thp = int(sys.argv[3]) if len(sys.argv)>3 else 1
a0 = 1/3
ch = build_chain(L, which_th=thp)
psi0 = W_state(L, a0)
F = factor_list(ch)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1', 'W2', 'W3AFF', 'W3AAF'))
W = ws.matrix()
np.save(f'W_{tag}.npy', W)
pickle.dump((F, ws.labels), open(f'W_{tag}.pkl', 'wb'))
print('saved', W.shape)
