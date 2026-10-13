import pickle, itertools
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
qs = pickle.load(open('qs.pkl', 'rb'))
F = {2: from_roots(thF), 3: from_roots(thF + [t + h for t in thF]),
     4: from_roots(thF + [t + h for t in thF] + [t + 2 * h for t in thF] + [t - h for t in thA])}
def Pcas(n, I):
    k = len(I)
    M = [[pscale(pshift(qs[n][m][0], -s * h), z[m - 1] ** (k - 1 - s)) for m in I] for s in range(k)]
    return ptrim(pdet(M))
allQ = []
maxrem = 0
for n in range(256):
    Q = {(i,): qs[n][i][0] for i in range(1, 5)}
    for k in (2, 3, 4):
        for I in itertools.combinations(range(1, 5), k):
            P = Pcas(n, I)
            q, r = pdivmod(P, F[k])
            rel = max([abs(x) for x in r] + [0]) / max(abs(x) for x in P)
            maxrem = max(maxrem, rel)
            Q[I] = ptrim(q)
    allQ.append(Q)
print('max relative remainder of Casoratians / vacuum factors', mp.nstr(maxrem, 5))
q1234 = [allQ[n][(1, 2, 3, 4)] for n in range(256)]
print('q1234 degrees', set(len(q) - 1 for q in q1234), 'values', mp.nstr(min(abs(q[0]) for q in q1234), 8), mp.nstr(max(abs(q[0]) for q in q1234), 8))
# degree table vs weight
from collections import Counter
c = Counter()
for n in range(256):
    w = states[n].weight
    degs = tuple(len(allQ[n][I]) - 1 for I in [(4,), (3, 4), (2, 3, 4), (1,), (1, 2), (1, 2, 3)])
    c[(w, degs)] += 1
for k in sorted(c)[:12]:
    print(k, c[k])
pickle.dump(allQ, open('allQ.pkl', 'wb'))
