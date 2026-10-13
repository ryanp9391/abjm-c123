import sympy as sp, pickle
th, h, R = sp.symbols('theta h Rhat')
cc = pickle.load(open('step1c.pkl', 'rb'))['coef']
v = sp.symbols('v1:5'); u = sp.symbols('u1:5')
for b in range(4):
    wb = sp.expand((-1)**b * sum(cc[(b, k)] * (h * th * R)**k for k in range(4)))
    # organise by v_i
    print('w_%d =' % b)
    for i in range(4):
        part = sp.factor(sp.expand(wb).coeff(v[i]))
        print('   v%d * [' % (i + 1), part, ']')
