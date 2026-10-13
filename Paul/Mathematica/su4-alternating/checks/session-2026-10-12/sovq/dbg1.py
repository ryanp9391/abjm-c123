import pickle
import mpmath as mp
from qfun import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
qs = pickle.load(open('qs.pkl', 'rb'))
for n, st in enumerate(states):
    if st.weight == (2, 0, 0, -2):
        print('vacuum index', n, [len(qs[n][i][0]) - 1 for i in range(1, 5)])
        u = mp.mpf('0.377')
        for i in range(1, 5):
            r = baxter_row(st, z[i - 1], u, 0)
            print(i, mp.nstr(r[0], 8))
        print([mp.nstr(st.tau(a, u), 10) for a in range(5)])
