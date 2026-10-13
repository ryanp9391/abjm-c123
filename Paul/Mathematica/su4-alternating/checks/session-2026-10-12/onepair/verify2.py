from common import *
thS, hS = sp.symbols('theta h')
cc = pickle.load(open('step1c.pkl', 'rb'))['coef']
vS = sp.symbols('v1:5'); uS = sp.symbols('u1:5')
cf = {k: sp.lambdify((thS, hS) + vS + uS, e, 'mpmath') for k, e in cc.items()}
random.seed(1)
maxdev = 0
for trial in range(3):
    v = [mp.mpf(random.randint(-9, 9)) / random.randint(1, 7) for _ in range(4)]
    u = [mp.mpf(random.randint(-9, 9)) / random.randint(1, 7) for _ in range(4)]
    c = {k: f(th, h, *v, *u) for k, f in cf.items()}
    for s in states:
        psi = s['ketm']
        lhs = sum(v[i] * u[j] * psi[4 * i + j] for i in range(4) for j in range(4)) / psi[3]
        r = tau(s, 1, th2)
        w = [(-1)**b * sum(c[(b, k)] * r**k for k in range(4)) for b in range(4)]
        Q123 = det([[QA(s, i, b) for b in (1, 2, 3)] for i in range(3)])
        rhs = det([w] + [[QA(s, i, b) for b in range(4)] for i in range(3)]) / Q123
        maxdev = max(maxdev, abs(lhs / rhs - 1))
print('det formula, 3 random (v,u) x 16 states: max rel dev', mp.nstr(maxdev, 3))
