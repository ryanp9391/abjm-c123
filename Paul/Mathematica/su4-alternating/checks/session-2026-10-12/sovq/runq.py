import pickle, time
import mpmath as mp
from qfun import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
t0 = time.time()
out = []
for n, st in enumerate(states):
    qs = {}
    for i in range(1, 5):
        q, r = solve_q(st, z[i - 1])
        qs[i] = (q, r)
    out.append(qs)
    if n % 32 == 0:
        print(n, st.weight, [len(qs[i][0]) - 1 for i in range(1, 5)], ['%.1e' % float(qs[i][1]) for i in range(1, 5)], '%.1f' % (time.time() - t0), flush=True)
pickle.dump(out, open('qs.pkl', 'wb'))
from collections import Counter
print(Counter(tuple(len(o[i][0]) - 1 for i in range(1, 5)) for o in out))
print('max rel sing', max(float(o[i][1]) for o in out for i in range(1, 5)))
