import pickle, itertools
import mpmath as mp
from collections import Counter
from qfun import *
from poly import *
mp.mp.dps = 50
data = pickle.load(open('eigen.pkl', 'rb')); states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
vac = [n for n, s in enumerate(states) if s.weight == (2, 0, 0, -2)][0]
print('vacuum index', vac)
for I in [(1,), (2,), (3,), (4,), (3, 4), (2, 3, 4), (1, 2, 3)]:
    print('vacuum q%s =' % (I,), [mp.nstr(c, 10) for c in allQ[vac][I]])
print('q1234 =', mp.nstr(allQ[0][(1, 2, 3, 4)][0], 20))
# degree table: (deg q4, q34, q234 ; q1, q12, q123) vs weight
c = Counter()
for n, s in enumerate(states):
    c[tuple(len(allQ[n][I]) - 1 for I in [(1,), (2,), (3,), (4,)])] += 1
print('deg (q1..q4) distribution:', len(c), 'patterns; all = 2 + mu:', all(tuple(len(allQ[n][(i,)]) - 1 for i in (1, 2, 3, 4)) == tuple(2 + m for m in s.weight) for n, s in enumerate(states)))
# nested degrees: deg q_{(4)}, q_{(34)}, q_{(234)} = K's ?
ok = True
for n, s in enumerate(states):
    mu = s.weight
    K = (2 + mu[3], 2 + mu[3] + 2 + mu[2] - 2, None)
print('example state 100 weight', states[100].weight, {I: len(allQ[100][I]) - 1 for I in allQ[100]})
# simple-spectrum statistics
ev = [s.d['combEig'] for s in states]
print('min gap', mp.nstr(min(abs(ev[i] - ev[j]) for i in range(256) for j in range(i)), 5))
print('# complex states', sum(1 for s in states if any(abs(mp.im(x)) > 1e-30 for x in s.tb[1])))
