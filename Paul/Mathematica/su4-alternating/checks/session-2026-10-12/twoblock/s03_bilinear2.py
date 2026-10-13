import sympy as sp
exec(open('s02_bilinear.py').read().split("print('=== (a) F-only")[0])
def show(name, K):
    d = sp.factor(sp.expand(K.det()))
    print(name); print('   K =', [[sp.factor(e) for e in row] for row in K.tolist()]); print('   det =', d)
# single-site cheap directions: site1 {e1,e4}, site2 {e4,e3}, site3 {e1,e2}, site4 {e4,e3}
show('(b) A-only: site 1 span(e1,e4) x site 3 span(e1,e2)', form(0, [0, 3], 2, [0, 1]))
show('(c) A3 span(e1,e2) x F2 span(e4,e3)', form(2, [0, 1], 1, [3, 2]))
show('(c) A3 span(e1,e2) x F4 span(e4,e3)', form(2, [0, 1], 3, [3, 2]))
show('(c) A1 span(e1,e4) x F2 span(e4,e3)', form(0, [0, 3], 1, [3, 2]))
show('(c) A1 span(e1,e4) x F4 span(e4,e3)', form(0, [0, 3], 3, [3, 2]))
# single-site expansions
for s, comps in [(0, [3, 1, 2]), (1, [2, 1, 0]), (2, [1, 2, 3]), (3, [2, 1, 0])]:
    for cc in comps:
        c = list(vac); c[s] = cc
        print(f' single site {s+1} component e{cc+1}: K =', sp.factor(Kpoly(tuple(c))))
