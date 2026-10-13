import pickle, itertools
import numpy as np
import mpmath as mp
from qfun import *
from poly import *
data = pickle.load(open('eigen.pkl', 'rb'))
states = [State(d) for d in data['states']]
allQ = pickle.load(open('allQ.pkl', 'rb'))
NSt = 256
def fit(target, cands, names):
    """least squares in log|.| : log|T| = c + sum m_j log|f_j|"""
    Y = np.array([float(mp.log(abs(t))) for t in target])
    X = np.array([[1.0] + [float(mp.log(abs(f[n]))) for f in cands] for n in range(NSt)])
    sol, res, rk, sv = np.linalg.lstsq(X, Y, rcond=None)
    resid = np.abs(X @ sol - Y).max()
    return sol, resid
QI = [(4,), (3, 4), (2, 3, 4), (1,), (1, 2), (1, 2, 3), (2,), (3,), (1, 3), (1, 4), (2, 3), (2, 4), (1, 2, 4), (1, 3, 4)]
def qval(X, u):
    return [peval(allQ[n][X], u) for n in range(NSt)]
# targets
tgt = {}
for s, t in [('F2', th[1]), ('F4', th[3])]:
    tgt[s + ' tau1(thF)'] = [st.tau(1, t) for st in states]
for s, t in [('A1', th[0]), ('A3', th[2])]:
    for a in (1, 2, 3):
        tgt[s + ' tau%d(thA-h)' % a] = [st.tau(a, t - h) for st in states]
shifts = range(-4, 4)
for name, T in tgt.items():
    site = name[:2]
    t = {'F2': th[1], 'F4': th[3], 'A1': th[0], 'A3': th[2]}[site]
    best = []
    for X in QI:
        cands = [qval(X, t + j * h) for j in shifts]
        sol, resid = fit(T, cands, None)
        best.append((resid, X, np.round(sol[1:], 3)))
    best.sort(key=lambda b: b[0])
    print(name)
    for b in best[:16]:
        print('   resid %.2e  Q%s  exps %s' % (b[0], b[1], b[2] if b[0] < 1e-8 else ''))
