import sympy as sp, pickle
th1, th2, h = sp.symbols('theta1 theta2 h')
d = pickle.load(open('step1.pkl', 'rb')); W = d['W']
lab = [(a, k) for a in range(4) for k in range(4)]
pb = [(i, j) for i in range(4) for j in range(4)]   # 0-based product basis e_{i+1} (x) e_{j+1}
for r, (a, k) in enumerate(lab):
    supp = [(i+1, j+1, sp.factor(W[r, 4*i+j])) for (i, j) in pb if W[r, 4*i+j] != 0]
    grades = set(i - j + 3 for (i, j, _) in supp)
    print((a, k), 'grade(s)', grades, supp)
