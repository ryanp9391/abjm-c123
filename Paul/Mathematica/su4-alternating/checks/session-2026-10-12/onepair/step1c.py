import sympy as sp, pickle
th1, th2, h, th = sp.symbols('theta1 theta2 h theta')
d = pickle.load(open('step1.pkl', 'rb')); W = d['W'].subs(th1, th2 + th)
W = W.applyfunc(sp.factor)
Winv = W.inv().applyfunc(sp.factor)
lab = [(a, k) for a in range(4) for k in range(4)]
v = sp.symbols('v1:5'); u = sp.symbols('u1:5')
vu = sp.Matrix([[v[i]*u[j] for i in range(4) for j in range(4)]])
cvec = (vu * Winv).applyfunc(sp.factor)
coef = {}
for r, (a, k) in enumerate(lab):
    coef[(a, k)] = sp.factor(cvec[0, r])
    print((a, k), ':', coef[(a, k)])
pickle.dump({'W': W, 'Winv': Winv, 'coef': coef}, open('step1c.pkl', 'wb'))
# check reconstruction
rec = sum((coef[lab[r]] * W[r, :] for r in range(16)), sp.zeros(1, 16))
print('reconstruction ok:', (rec - vu).applyfunc(sp.simplify) == sp.zeros(1, 16))
