import sys, pickle
from pipeline import *
from span import *
from sparse import *
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
a0 = 1/3
ch = build_chain(L)
psi0 = W_state(L, a0)
F = factor_list(ch)
n = len(F)
ws, B1, B2 = build_positive_words(ch, psi0, F, ('W1', 'W2'), verbose=False)
Wpos = ws.matrix(); labs_pos = ws.labels
nsing = 1 + n
# inverse words: all invertible inverses x base words with <= 2 positive factors? for L=2 do base = id + singles
base = Wpos[:, :nsing]; base_l = labs_pos[:nsing]
Vinv1, Linv1, info = build_inverse_words(ch, F, range(n), base, base_l, verbose=False)
print('inverse condition (pivot ratios) min', min(i[1] for i in info))
w0 = W_state_rest(L, a0)
T = targets(L)
def report(name, Wm):
    Q, gaps = orth_basis(Wm)
    s, U, Vt = localised(Q, w0, None)
    nloc = int((s > 1 - 1e-8).sum())
    out = [name, Wm.shape[1], Q.shape[1], nloc]
    out += ['%.1e' % resid(Q, tv) for tv in T.values()]
    return out, Q

named = [factor_index(F, 1, 0, 0), factor_index(F, 1, 1, 0)]   # t1(th1), t1(th2)
nB = nsing
def inv_cols(idx_list, base_range):
    # Vinv1 blocks are in order of i in range(n), each of nB columns (skipped near-singular ones would shift; none skipped at these params)
    cols = []; 
    for i in idx_list:
        for b in base_range:
            cols.append(i * nB + b)
    return Vinv1[:, cols]
sets = {}
sets['S1: 1+singles'] = Wpos[:, :nsing]
sets['S1i2: S1 + t1(th1)^-1,t1(th2)^-1'] = np.concatenate([Wpos[:, :nsing], inv_cols(named, [0])], axis=1)
sets['S1iAll: S1 + all inverse singles'] = np.concatenate([Wpos[:, :nsing], inv_cols(range(n), [0])], axis=1)
sets['S2: 1+singles+pairs'] = Wpos
sets['S2i2: S2 + t1(th1)^-1, t1(th2)^-1 x (1,singles)'] = np.concatenate([Wpos, inv_cols(named, range(nB))], axis=1)
sets['S2iAll: S2 + all inverse x (1,singles)'] = np.concatenate([Wpos, inv_cols(range(n), range(nB))], axis=1)
rows = []
for k, Wm in sets.items():
    r, Q = report(k, Wm); rows.append(r)
print(['name','Nw','rank','nloc'] + list(T.keys()))
for r in rows: print(r)
