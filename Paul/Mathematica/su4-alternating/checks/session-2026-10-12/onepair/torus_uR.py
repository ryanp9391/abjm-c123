import sympy as sp
cs = sp.symbols('c1:5'); k, R, y = sp.symbols('kappa R y')
G = sp.Matrix([[cs[0], -cs[1], cs[2], -cs[3]], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
chi = lambda x: x**4 - cs[0]*x**3 + cs[1]*x**2 - cs[2]*x + cs[3]
u = sp.Matrix([[0, 0, 0, 1]]) * (sp.eye(4) - k * G).inv()
uR = sp.simplify((u * sp.Matrix([R**3, R**2, R, 1]))[0])
print('u(R) == (1 - chi(R)/chi(1/kappa))/(1 - kappa R):', sp.simplify(uR - (1 - chi(R) / chi(1 / k)) / (1 - k * R)) == 0)
v = sp.Matrix([[1, 0, 0, 0]]) * (sp.eye(4) - k * G.T)
print('v = e1(1 - kappa G^T) =', list(v))
