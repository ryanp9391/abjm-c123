import sympy as sp, pickle
th, h, R = sp.symbols('theta h R')
cc = pickle.load(open('step1c.pkl', 'rb'))['coef']
v = sp.symbols('v1:5'); u = sp.symbols('u1:5')
uR = sum(u[j] * R**(3 - j) for j in range(4))
S = sum(v[i] * u[j] * R**(i - j) for i in range(4) for j in range(4) if i >= j)
T = sum(v[i] * u[j] * R**(3 + i - j) for i in range(4) for j in range(4) if i < j)
N = [sp.factorial(b) * h**b * sp.prod([(m * h - th) for m in range(1, b + 2)]) for b in range(4)]
ok1 = ok2 = True
for b in range(4):
    wb = (-1)**b * sum(cc[(b, k)] * (h * th * R)**k for k in range(4))
    nu = (b * h - th) * v[b] + h * sum(v[i] * R**(i - b) for i in range(b, 4))
    num1 = uR * nu - h * R**(3 - b) * S
    ok1 &= sp.simplify(wb - num1 / N[b]) == 0
    nup = (b * h - th) * v[b] - h * sum(v[i] * R**(i - b) for i in range(0, b))   # i (0-based) <= b-1  -> R^{i-1-b} with 1-based
    num2 = uR * nup + h * R**(-b) * T
    ok2 &= sp.simplify(wb - num2 / N[b]) == 0
print('closed form 1 (S, includes v.u):', ok1)
print('closed form 2 (T, strictly lower):', ok2)
