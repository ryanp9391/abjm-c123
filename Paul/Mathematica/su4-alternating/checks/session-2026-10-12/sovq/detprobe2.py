import pickle, itertools, sys
import numpy as np
import mpmath as mp
from sov_l4 import WORDS, TH0, H0
Mt = pickle.load(open('Mt.pkl', 'rb'))
th = [float(t) for t in TH0]; h = float(H0)
keys = list(Mt.keys()); vals = np.array([complex(Mt[k]) for k in keys])
yi = np.array([WORDS.index(k[0]) for k in keys]); xi = np.array([WORDS.index(k[1]) for k in keys])
rng = np.random.default_rng(1)
def rand_data():
    c = lambda *s: rng.normal(size=s) + 1j * rng.normal(size=s)
    fF = {s: c(2) for s in (1, 3)}            # f(th-h), f(th)
    fA = {s: c(3, 5) for s in (0, 2)}         # f_i(th - j h), j = 0..4 (use 1..4)
    gA = {s: c(2) for s in (0, 2)}            # g(th-h), g(th-2h)
    gF = {s: c(3, 4) for s in (1, 3)}         # g_i(th + (j-1)h), j = 0..3
    return fF, fA, gA, gF
def B(d):
    fF, fA, gA, gF = d
    Fv = np.ones(256, complex); Gv = np.ones(256, complex)
    for n, x in enumerate(WORDS):
        for s in (0, 2):
            cols = [j for j in (1, 2, 3, 4) if j != x[s] + 1]
            Fv[n] *= np.linalg.det(fA[s][:, cols])
        for s in (1, 3):
            Fv[n] *= fF[s][0] ** (3 - x[s]) * fF[s][1] ** x[s]
        y = x
        for s in (0, 2):
            Gv[n] *= gA[s][0] ** (3 - y[s]) * gA[s][1] ** y[s]
        for s in (1, 3):
            cols = [j for j in (0, 1, 2, 3) if j != y[s]]
            Gv[n] *= np.linalg.det(gF[s][:, cols])
    return np.sum(Gv[yi] * vals * Fv[xi])
def rowparts(d, variant):
    """for each row (a,s): two 12-vectors (p0 part, p1 part); columns (b, beta)"""
    fF, fA, gA, gF = d
    rows = []
    for s in range(4):
        for a in range(3):
            parts = []
            for p in (0, 1):
                vec = np.zeros(12, complex)
                for b in range(3):
                    for beta in range(4):
                        if s in (1, 3):
                            P = th[s] - h + p * h                     # bra points th-h, th
                            jg = p + b if variant[0] == '+' else p + 2 - b
                            val = fF[s][p] * gF[s][a, jg]
                        else:
                            P = th[s] - h - p * h                     # ket points th-h, th-2h
                            jf = 1 + p + b if variant[1] == '+' else 1 + p + 2 - b
                            val = gA[s][p] * fA[s][a, jf]
                        vec[3 * beta + b] = P ** beta * val
                parts.append(vec)
            rows.append((s, parts))
    return rows
def coeffs(d, variant):
    """c_k for k in {0..3}^4 : det E = sum_k c_k rho^k"""
    rows = rowparts(d, variant)
    c = np.zeros((4, 4, 4, 4), complex)
    choices = np.array(list(itertools.product((0, 1), repeat=12)))
    mats = np.empty((len(choices), 12, 12), complex)
    for r, (s, parts) in enumerate(rows):
        mats[:, r, :] = np.where(choices[:, r:r + 1] == 0, parts[0][None, :], parts[1][None, :])
    dets = np.linalg.det(mats)
    for ch, dv in zip(choices, dets):
        k = tuple(int(ch[3 * s:3 * s + 3].sum()) for s in range(4))
        c[k] += dv
    return c.reshape(-1)
for variant in ['++', '+-', '-+', '--']:
    ndata = 300
    A = []; bb = []
    for _ in range(ndata):
        d = rand_data()
        A.append(coeffs(d, variant)); bb.append(B(d))
    A = np.array(A); bb = np.array(bb)
    sol, res, rk, sv = np.linalg.lstsq(A, bb, rcond=None)
    rel = np.linalg.norm(A @ sol - bb) / np.linalg.norm(bb)
    print('variant', variant, 'rank', rk, 'relative residual of linear fit', rel)
    u = sol.reshape(4, 4, 4, 4)
    # product-form test: u_k = C prod rho_s^{k_s}
    C = u[0, 0, 0, 0]
    rho = [u[tuple(1 if t == s else 0 for t in range(4))] / C for s in range(4)]
    pred = np.array([[[[C * rho[0] ** a * rho[1] ** b * rho[2] ** c * rho[3] ** e for e in range(4)] for c in range(4)] for b in range(4)] for a in range(4)])
    print('   C =', C, ' rho =', np.round(rho, 8), ' product-form defect', np.abs(pred - u).max() / np.abs(u).max())
    sys.stdout.flush()
