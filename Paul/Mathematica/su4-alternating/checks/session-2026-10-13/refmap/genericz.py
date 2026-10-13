"""genericz.py -- K_kappa admissibility at GENERIC z = (3/2, 5/4, 2/3, 4/7) (z2 != z3; ABJM-point words, pair level): enumerate the zero patterns
(supports of v and u) with generic values on the support and v.u = 0 imposed, test T_kappa in the 5-word span (exact, 3 kappas).  Output: genericz_out.txt"""
from helpers import *
import itertools, random
out = open('genericz_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
zg = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
random.seed(7)
def rnd(): return Fr(random.randint(1, 9) * random.choice([-1, 1]), random.randint(1, 5))
def sample(Sv, Su):
    v = [rnd() if j in Sv else Fr(0) for j in range(4)]; u = [rnd() if j in Su else Fr(0) for j in range(4)]
    common = [j for j in range(4) if j in Sv and j in Su]
    if common:
        if len(common) == 1:
            return None   # v.u = v_j u_j != 0 forced
        j = common[0]; u[j] = u[j] - sum(a * b for a, b in zip(v, u)) / v[j]
        if u[j] == 0: return None
    else:
        pass   # disjoint supports: v.u = 0 automatically
    return v, u
res = {}
for rv in range(1, 5):
    for Sv in itertools.combinations(range(4), rv):
        for ru in range(1, 5):
            for Su in itertools.combinations(range(4), ru):
                s = sample(set(Sv), set(Su))
                if s is None: continue
                v, u = s
                if sum(a * b for a, b in zip(v, u)) != 0: continue
                ok = test(v, u, zg)[2]
                res[(Sv, Su)] = ok
adm = [k for k, ok in res.items() if ok]
P('generic z: patterns tested', len(res), ' admissible', len(adm))
nm = lambda S: ''.join(str(j + 1) for j in S)
P('admissible (supp v ; supp u):', ['%s;%s' % (nm(a), nm(b)) for a, b in adm])
P('NOT admissible patterns with |supp v|+|supp u| <= 5:', ['%s;%s' % (nm(a), nm(b)) for (a, b), ok in res.items() if not ok and len(a) + len(b) <= 5])
out.close()
