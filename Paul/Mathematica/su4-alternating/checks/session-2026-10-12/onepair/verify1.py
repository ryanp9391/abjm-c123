from common import *
maxdev = {'r': 0, 'tauA123': 0, 'tauA_other': 0, 'tau4_A': 0, 'tF_ge2': 0, 'Q1234(th1)': 0}
for s in states:
    r = tau(s, 1, th2)
    for j in range(4):
        maxdev['r'] = max(maxdev['r'], abs(QF(s, j, 1) / QF(s, j, 0) / r - 1))
    for a in (2, 3, 4): maxdev['tF_ge2'] = max(maxdev['tF_ge2'], abs(tau(s, a, th2)))
    maxdev['tau4_A'] = max(maxdev['tau4_A'], abs(tau(s, 4, th1 - h)))
    for trip in combinations(range(4), 3):
        D = [det([[QA(s, i, b) for b in range(4) if b != a] for i in trip]) for a in range(4)]
        for a in range(4):
            dev = abs(D[a] / D[0] / tau(s, a, th1 - h) - 1)
            key = 'tauA123' if trip == (0, 1, 2) else 'tauA_other'
            maxdev[key] = max(maxdev[key], dev)
    q4 = det([[QA(s, i, b) for b in range(4)] for i in range(4)])
    maxdev['Q1234(th1)'] = max(maxdev['Q1234(th1)'], abs(q4) / mp.fprod([abs(QA(s, i, 0)) for i in range(4)]))
for k, v in maxdev.items(): print(k, mp.nstr(v, 3))
