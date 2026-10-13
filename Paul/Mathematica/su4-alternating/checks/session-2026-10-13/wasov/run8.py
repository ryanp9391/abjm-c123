from wa_sov_core import *
cn = pickle.load(open(os.path.join(CACHE, 'cn.pkl'), 'rb'))
idx = {w: i for i, w in enumerate(WORDS)}
for n in range(5):
    m = 4 - n
    for name, mapx in (('x -> 3-x', lambda x: tuple(3 - c for c in x)), ('identity', lambda x: x), ('A:3-a, F:k', lambda x: (3 - x[0], x[1], 3 - x[2], x[3])), ('A:a, F:3-k', lambda x: (x[0], 3 - x[1], x[2], 3 - x[3]))):
        ratio = [cn[m][idx[mapx(x)]] / cn[n][idx[x]] for x in WORDS]
        print('c^(%d)_{map(x)}/c^(%d)_x with map = %-12s: #distinct ratios %d, tensor ranks %s' % (m, n, name, len(set(ratio)), tensor_ranks(ratio)['site']))
