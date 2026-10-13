import sympy as sp, mpmath as mp, pickle
exec(open('L2num2.py').read().split('results = []')[0])
d0 = 2*h; u0 = d0 - h; thA, thF = d0, 0
zs = [mp.mpf(3)/2, mp.mpf(5)/4, mp.mpf(2)/3, mp.mpf(4)/7]
# polynomial coefficient matrices of t_a(uu) at d = d0
Cm = {}
for a in range(1,5):
    P = sp.Poly(tsym[a].subs(d, d0).applyfunc(sp.expand)[0,0]*0 + 0, uu)  # dummy
    mats = []
    deg = 2*a
    Ma = tsym[a].subs(d, d0)
    for m in range(deg+1):
        mats.append(tomp(Ma.applyfunc(lambda e: sp.Poly(e, uu).coeff_monomial(uu**m) if e != 0 else 0)))
    Cm[a] = mats
def taupoly(a, psi):
    j = max(range(16), key=lambda i: abs(psi[i]))
    return [ (M*psi)[j]/psi[j] for M in Cm[a] ]
def pev(c, x): return mp.fsum(cc*x**m for m, cc in enumerate(c))
def pder(c, x): return mp.fsum(m*cc*x**(m-1) for m, cc in enumerate(c) if m)
V4 = lambda x: (x - thA)*(x - thF - h)
def Pa(a, x):
    r = mp.mpf(1)
    for bb in range(a, 4): r *= V4(x - bb*h)
    return r
def Pa_d(a, x):
    return mp.diff(lambda y: Pa(a, y), x)
def solve_q(taus, zeta, Dmax=8):
    # find polynomial q of minimal degree with B(u) = 0
    for D in range(Dmax+1):
        # sample B at many points -> linear system in q coefficients
        npts = 8 + D + 6
        pts = [mp.mpf(k)/7 + mp.mpf(1)/13 for k in range(npts)]
        A = mp.matrix(npts, D+1)
        for r, x in enumerate(pts):
            for m in range(D+1):
                val = 0
                for a in range(5):
                    ta = 1 if a == 0 else pev(taus[a], x)
                    val += (-1)**a * ta * zeta**(1-a) * (x + h - a*h)**m * Pa(a, x)
                A[r, m] = val
        U, S, V = mp.svd_r(A) if False else (None, None, None)
        # nullspace via QR-like: use mp.svd
        U, S, V = mp.svd(A)
        smin = S[D] if len(S) > D else 0
        if abs(smin) < mp.mpf(10)**(-40) * max(1, abs(S[0])):
            q = [V[D, m] for m in range(D+1)]
            return q, D
    return None, None
psis = eigstates(d0)
report = []
for n, psi in enumerate(psis):
    taus = {a: taupoly(a, psi) for a in range(1,5)}
    # direct eigenvalues
    x_dir = pev(taus[1], u0); tF = pev(taus[1], thF)
    w_dir = pder(taus[2], u0) - pder(taus[1], u0)*tF; y_dir = pder(taus[3], u0)
    chk_t4 = abs(pev(taus[4], u0)); chk_t2F = max(abs(pev(taus[a], thF)) for a in (2,3,4))
    best = None
    for conv in ('z', '1/z'):
        qs = []
        for zj in zs:
            zeta = zj if conv == 'z' else 1/zj
            q, D = solve_q(taus, zeta)
            qs.append((q, D, zeta))
        if all(q is not None for q,_,_ in qs): best = (conv, qs); break
    conv, qs = best
    # F relation
    Frel = max(abs(zeta*V4(thF)*pev(q, thF+h)/pev(q, thF) - tF) for q,_,zeta in qs)
    # derivative Baxter: rows j
    rows = []; rhs = []
    for q, D, zeta in qs:
        qf = lambda x: pev(q, x)
        Rj = -mp.diff(lambda y: pev(q, y)*Pa(1, y), u0) + tF/zeta*mp.diff(lambda y: pev(q, y-h)*Pa(2, y), u0)
        rows.append([Rj, Pa(2,u0)*qf(thF)/zeta, -Pa(3,u0)*qf(thF-h)/zeta**2])
        rhs.append(-zeta*qf(thF+2*h)*Pa_d(0,u0) - pder(taus[4], u0)*qf(thF-2*h)/zeta**3)
    A3 = mp.matrix(rows[:3]); b3 = mp.matrix(rhs[:3])
    sol = mp.lu_solve(A3, b3)
    resid4 = abs(mp.fsum(rows[3][i]*sol[i] for i in range(3)) - rhs[3])
    errs = (abs(sol[0]-x_dir), abs(sol[1]-w_dir), abs(sol[2]-y_dir))
    report.append((n, conv, [D for _,D,_ in qs], mp.nstr(chk_t4,3), mp.nstr(chk_t2F,3), mp.nstr(Frel,3), [mp.nstr(e,3) for e in errs], mp.nstr(resid4,3), mp.nstr(x_dir,8), mp.nstr(w_dir,8), mp.nstr(y_dir,8)))
for r in report: print(r)
