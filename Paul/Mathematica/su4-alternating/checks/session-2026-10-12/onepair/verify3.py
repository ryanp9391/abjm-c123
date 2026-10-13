# polynomial-twisted gauge form, decomposition (product channel + R-Casoratian channel), chi_A(R)
from common import *
def Qt(s, j, x): return z[j]**(x / h) * qj(s, j, x)   # Q-tilde = z^{u/h} q(u)
xs = [th1 - b * h for b in range(4)]
random.seed(7)
dev = {'Qt-form': 0, 'decomp': 0, 'R=Qt ratio': 0}
casR = {}
for trial in range(3):
    v = [mp.mpf(random.randint(-9, 9)) / random.randint(1, 7) for _ in range(4)]
    u = [mp.mpf(random.randint(-9, 9)) / random.randint(1, 7) for _ in range(4)]
    for s in states:
        psi = s['ketm']
        lhs = sum(v[i] * u[j] * psi[4 * i + j] for i in range(4) for j in range(4)) / psi[3]
        R = tau(s, 1, th2) / (h * th)
        for j in range(4): dev['R=Qt ratio'] = max(dev['R=Qt ratio'], abs(Qt(s, j, th2 + h) / Qt(s, j, th2) / R - 1))
        uR = sum(u[j] * R**(3 - j) for j in range(4))
        S = sum(v[i] * u[j] * R**(i - j) for i in range(4) for j in range(4) if i >= j)
        T = sum(v[i] * u[j] * R**(3 + i - j) for i in range(4) for j in range(4) if i < j)
        num = [uR * ((b * h - th) * v[b] + h * sum(v[i] * R**(i - b) for i in range(b, 4))) - h * R**(3 - b) * S for b in range(4)]
        Qrows = [[Qt(s, i, x) for x in xs] for i in range(3)]
        Q123 = det([row[1:] for row in Qrows])
        rhs = det([num] + Qrows) / ((h - th) * Q123)
        dev['Qt-form'] = max(dev['Qt-form'], abs(lhs / rhs - 1))
        nup = [(b * h - th) * v[b] - h * sum(v[i] * R**(i - b) for i in range(0, b)) for b in range(4)]
        cas = det([[R**(-b) for b in range(4)]] + Qrows)
        rhs2 = (uR * det([nup] + Qrows) + h * T * cas) / ((h - th) * Q123)
        dev['decomp'] = max(dev['decomp'], abs(lhs / rhs2 - 1))
        casR[s['label']] = cas / Q123
for k, x in dev.items(): print(k, mp.nstr(x, 3))
print('R-Casoratian det[R^{-b}; Qt1..3]/Qt123(th1-h) per state (= R^-3 chi_A(R)):')
for k, x in casR.items(): print('  ', k, mp.nstr(x, 6))
