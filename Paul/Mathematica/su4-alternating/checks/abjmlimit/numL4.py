"""High-precision (mpmath) check at L = 4 on twisted eigenstates, from exact rational matrices (crt builds)."""
import pickle, sys, time
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
from words import WIDX

mp.mp.dps = 60


def to_mp(Mx):
    return mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in Mx])


def to_np(Mx):
    return np.array([[float(x) for x in row] for row in Mx])


def bra_exact(vs):
    out = [Fr(1)]
    for v in vs:
        out = [a * Fr(b) for a in out for b in v]
    return out


B_TEST = bra_exact([[1, Fr(2, 3), Fr(-1, 2), Fr(3, 5)], [Fr(2, 7), -1, Fr(1, 3), 1], [Fr(-4, 3), 1, Fr(5, 2), Fr(1, 9)], [1, Fr(3, 4), Fr(-2, 5), Fr(7, 3)]])


def refine(Amp, Anp, lam0, x0, iters=10):
    n = Anp.shape[0]
    pidx = int(np.argmax(np.abs(x0)))
    x0 = x0 / x0[pidx]
    x = mp.matrix([mp.mpc(complex(v)) for v in x0])
    lam = mp.mpc(complex(lam0))
    J = np.zeros((n + 1, n + 1), dtype=complex)
    J[:n, :n] = Anp - lam0 * np.eye(n)
    J[:n, n] = -x0
    J[n, pidx] = 1
    lu = np.linalg.inv(J)
    res = None
    for it in range(iters):
        r = Amp * x - lam * x
        res = max(abs(v) for v in r)
        if res < mp.mpf(10) ** (-(mp.mp.dps - 8)):
            break
        rn = np.array([complex(v) for v in r] + [0j])
        dz = -lu @ rn
        for i in range(n):
            x[i] += mp.mpc(dz[i])
        lam += mp.mpc(dz[n])
    return lam, x, res


def mvec(Mmp, x):
    return Mmp * x


def run(tag, nstates=8):
    t0 = time.time()
    M = pickle.load(open('exact_M_%s.pkl' % tag, 'rb'))
    T = pickle.load(open('exact_t1_%s.pkl' % tag, 'rb'))
    Mmp = to_mp(M)
    bmp = mp.matrix([mp.mpf(x.numerator) / x.denominator for x in B_TEST])
    c = mp.lu_solve(Mmp.T, bmp)
    t1 = time.time() - t0
    Tmp = to_mp(T)
    Tnp = to_np(T)
    ev, evec = np.linalg.eig(Tnp)
    # choose well separated eigenvalues
    gaps = np.array([np.min(np.abs(np.delete(ev, i) - ev[i])) for i in range(len(ev))])
    order = np.argsort(-gaps)[:nstates]
    rows = []
    idx0 = WIDX.index((0, 0, 0, 0))
    for i in order:
        lam, psi, res = refine(Tmp, Tnp, ev[i], evec[:, i])
        ov = Mmp * psi
        o0 = ov[idx0]
        tA1 = [ov[WIDX.index((a, 0, 0, 0))] / o0 for a in range(4)]
        tA3 = [ov[WIDX.index((0, 0, a, 0))] / o0 for a in range(4)]
        tF2 = ov[WIDX.index((0, 1, 0, 0))] / o0
        tF4 = ov[WIDX.index((0, 0, 0, 1))] / o0
        pred = [tA1[w[0]] * tF2 ** w[1] * tA3[w[2]] * tF4 ** w[3] for w in WIDX]
        fact = max(abs(ov[j] / o0 - pred[j]) / max(1, abs(pred[j])) for j in range(256))
        terms = [c[j] * pred[j] for j in range(256)]
        S = mp.fsum(terms)
        direct = mp.fsum(bmp[j] * psi[j] for j in range(256)) / o0
        rows.append((mp.nstr(res, 3), mp.nstr(fact, 3), mp.nstr(S, 12), mp.nstr(direct, 12), mp.nstr(abs(S - direct) / abs(direct), 3),
                     mp.nstr(max(abs(t) for t in terms) / abs(S), 3)))
    maxc = max(abs(v) for v in c)
    out = {'tag': tag, 'maxc': maxc, 'c': c, 'rows': rows, 'time': time.time() - t0, 'solve_time': t1}
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    out = run(tag)
    print('delta = 1/%s  max|c| =' % tag, mp.nstr(out['maxc'], 6), ' solve', round(out['solve_time'], 1), 's, total', round(out['time'], 1), 's')
    for r in out['rows']:
        print('  eig-res %s  factorisation %s  sum c*tau = %s  direct = %s  rel.diff %s  max|term|/|sum| %s' % r)
    pickle.dump({'tag': tag, 'c': [complex(x) if False else str(x) for x in out['c']], 'rows': out['rows'], 'maxc': str(out['maxc'])}, open('numL4_%s.pkl' % tag, 'wb'))
