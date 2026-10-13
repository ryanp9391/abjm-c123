import sympy as sp, pickle
th, h = sp.symbols('theta h'); cs = sp.symbols('c1:5'); al, be, ga, de = sp.symbols('alpha beta gamma delta')
cc = pickle.load(open('step1c.pkl', 'rb'))['coef']
vS = sp.symbols('v1:5'); uS = sp.symbols('u1:5')
G = sp.Matrix([[cs[0], -cs[1], cs[2], -cs[3]], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
e1 = sp.Matrix([[1, 0, 0, 0]]); e4 = sp.Matrix([[0, 0, 0, 1]])
def supp(v, u):
    sub = dict(zip(vS, list(v))); sub.update(dict(zip(uS, list(u))))
    return {k: sp.factor(e.subs(sub)) for k, e in cc.items() if sp.simplify(e.subs(sub)) != 0}
cases = {
 'A only, deg 1: v=e1(al+be G^T), u=e4': (e1 * (al * sp.eye(4) + be * G.T), e4),
 'F only, deg 1: v=e1, u=e4(ga+de G)': (e1, e4 * (ga * sp.eye(4) + de * G)),
 'both deg 1': (e1 * (al * sp.eye(4) + be * G.T), e4 * (ga * sp.eye(4) + de * G)),
 'A only generic v, u=e4': (sp.Matrix([vS]), e4),
 'F only generic u, v=e1': (e1, sp.Matrix([uS])),
}
for name, (v, u) in cases.items():
    s = supp(v, u)
    print(name, ': #words', len(s))
    for k in sorted(s): print('    ', k, s[k])
# BPS: c30 simplification
v = sp.Matrix([vS]); u4b = -(vS[0]*uS[0] + vS[1]*uS[1] + vS[2]*uS[2]) / vS[3]
print('BPS c30 =', sp.factor(cc[(3, 0)].subs(uS[3], u4b)))
