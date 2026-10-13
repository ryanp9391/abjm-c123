# Q-functions (RV Gamma gauge) for all 16 states, numeric mpmath 50 digits
import sympy as sp, pickle, mpmath as mp
mp.mp.dps = 50
def M(x):
    if isinstance(x, sp.Basic):
        x = sp.nsimplify(x) if not x.is_Rational else x
        return mp.mpf(x.p) / x.q
    if isinstance(x, mp.mpc): return mp.re(x) if abs(mp.im(x)) < mp.mpf(10)**-40 else x
    return x
h = mp.mpf(1) / 2; th1 = mp.mpf(2) / 7; th2 = mp.mpf(-3) / 5
z = [mp.mpf(3) / 2, mp.mpf(5) / 4, mp.mpf(2) / 3, mp.mpf(4) / 7]
th = th1 - th2
V4 = lambda x: (x - th1) * (x - th2 - h)
states = pickle.load(open('states.pkl', 'rb'))
def taufun(s, a):
    if a == 0: return lambda x: mp.mpf(1)
    cs = [M(c) for c in s['tau'][a]]
    return lambda x: sum(cs[m] * x**m for m in range(len(cs)))
def solve_q(s, zj, D=3):
    taus = [taufun(s, a) for a in range(5)]
    pts = [mp.mpf(k) / 13 + mp.mpf(1) / 3 for k in range(-12, 12)]
    rows = []
    for x in pts:
        row = []
        for m in range(D + 1):
            row.append(sum((-1)**a * taus[a](x) * zj**(4 - a) * (x + h - a * h)**m * mp.fprod([V4(x - b * h) for b in range(a, 4)]) for a in range(5)))
        rows.append(row)
    A = mp.matrix(rows)
    Uu, S, Vt = mp.svd_r(A)
    vec = [Vt[D, m] for m in range(D + 1)]   # smallest singular value (svd_r returns descending)
    # normalise: leading nonzero coefficient = 1
    lead = max(m for m in range(D + 1) if abs(vec[m]) > mp.mpf(10)**-30)
    vec = [c / vec[lead] for c in vec]
    return vec, S[D], S[D - 1], lead
for s in states:
    s['tauf'] = [taufun(s, a) for a in range(5)]
    s['q'] = []
    info = []
    for j in range(4):
        vec, smin, snext, deg = solve_q(s, z[j])
        s['q'].append(vec)
        info.append((deg, mp.nstr(smin, 3), mp.nstr(snext, 3)))
    print(s['label'], info)
pickle.dump([{k: v for k, v in s.items() if k != 'tauf'} for s in states], open('qfun.pkl', 'wb'))
