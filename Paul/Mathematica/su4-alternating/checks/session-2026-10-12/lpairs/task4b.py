import itertools
import numpy as np, pickle
import lp, qdict
d = pickle.load(open('task4a.pkl', 'rb'))
h = d['h']; z = d['z']; ch = lp.chain_abjm(2, d['thF'], h, z, 'c')
worst = dict(tF=0, x=0, t2=0, t3=0)
for s in range(256):
    qs = d['Q'][s]
    for I in itertools.combinations(range(4), 3) if False else [(0, 1, 2), (1, 2, 3), (0, 1, 3), (0, 2, 3)]:
        tauF, x, t2, t3h, res = qdict.formulas(ch, qs, z, 1, I)
        rel = lambda a, b: abs(a - b) / max(abs(b), 1e-12)
        worst['tF'] = max(worst['tF'], max(rel(t, d['tF'][s]) for t in tauF))
        worst['x'] = max(worst['x'], rel(x, d['x'][s])); worst['t2'] = max(worst['t2'], rel(t2, d['t2u0'][s])); worst['t3'] = max(worst['t3'], rel(t3h, d['t3h'][s]))
    if s < 3: print(s, 'tauF', d['tF'][s], tauF[0], ' x', d['x'][s], x, ' t3h', d['t3h'][s], t3h, 'qI-fit res %.1e' % res)
print('max relative deviations over 256 states x 4 triples:', {k: '%.1e' % v for k, v in worst.items()})
