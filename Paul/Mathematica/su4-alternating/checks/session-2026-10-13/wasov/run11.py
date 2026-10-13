"""minimal-degree polynomial bases of the column/row space of the rank-2 matrices m^(n)(r2,r4), n = 0,1,3,4: m = U C V^T with C constant?"""
import pickle, itertools, os, sys, sympy as sp
WORDS = list(itertools.product(range(4), repeat=4))
cn = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'cn.pkl'), 'rb'))
r2, r4 = sp.symbols('r2 r4')
def mmat(c):
    c4 = {w: c[i] for i, w in enumerate(WORDS)}
    return sp.Matrix(4, 4, lambda a1, a3: sum(sp.Rational(c4[(a1, k2, a3, k4)].numerator, c4[(a1, k2, a3, k4)].denominator) * r2 ** k2 * r4 ** k4 for k2 in range(4) for k4 in range(4)))
def poly_kernel(m):
    """polynomial left-kernel vectors (n^T m = 0), cleared of denominators"""
    ns = m.T.nullspace()
    out = []
    for v in ns:
        den = sp.lcm([sp.fraction(sp.together(x))[1] for x in v])
        w = [sp.cancel(x * den) for x in v]
        g = sp.gcd(sp.gcd(w[0], w[1]), sp.gcd(w[2], w[3]))
        out.append([sp.cancel(x / g) for x in w])
    return out
def space_basis(kern, d, e):
    """polynomial vectors p of bidegree <= (d,e) with kern_k . p = 0 for all k"""
    mons = [(i, j) for i in range(d + 1) for j in range(e + 1)]
    cs = [[sp.Symbol('c_%d_%d_%d' % (a, i, j)) for (i, j) in mons] for a in range(4)]
    p = [sum(c * r2 ** i * r4 ** j for c, (i, j) in zip(cs[a], mons)) for a in range(4)]
    eqs = []
    for n in kern:
        expr = sp.expand(sum(n[a] * p[a] for a in range(4)))
        eqs += sp.Poly(expr, r2, r4).coeffs()
    unk = [c for row in cs for c in row]
    A = sp.Matrix([[sp.diff(eq, u) for u in unk] for eq in eqs])
    ns = A.nullspace()
    return [[sum(v[a * len(mons) + t] * r2 ** mons[t][0] * r4 ** mons[t][1] for t in range(len(mons))) for a in range(4)] for v in ns]
for n in (1, 0):
    m = mmat(cn[n])
    kl = poly_kernel(m)            # left kernel: n^T m = 0 => column space of m = {p: n.p = 0}
    kr = poly_kernel(m.T)          # right kernel => row space
    degl = [(max(sp.Poly(x, r2, r4).degree(r2) for x in k if x != 0), max(sp.Poly(x, r2, r4).degree(r4) for x in k if x != 0)) for k in kl]
    print('component', n, 'left-kernel bidegrees', degl, flush=True)
    for d in range(0, 4):
        Bu = space_basis(kl, d, d)
        print('  column space: # polynomial vectors of bidegree <= (%d,%d): %d' % (d, d, len(Bu)), flush=True)
        if len(Bu) >= 1 and d <= 2: pass
    # build U from degree-minimal vectors
    sys.stdout.flush()
