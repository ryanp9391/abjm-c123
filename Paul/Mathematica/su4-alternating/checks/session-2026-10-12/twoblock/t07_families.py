import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; rng = np.random.default_rng(5)
def rv(n): return [X.rnd(rng) for _ in range(n)]
def bps_pair_random():
    v = np.array(rv(4)); u = np.array(rv(4))
    # make v.u = 0 by adjusting u[3]
    s = int((v[:3] * u[:3] % p).sum() % p)
    u[3] = (-s * pow(int(v[3]), -1, p)) % p
    return v % p, u % p
def report(name, vs):
    b = X.bra(vs); c = X.coeffs(b)
    print(f'{name:70s}', X.ranks(c))
# companion-frame families
e1, e2, e3, e4 = [np.array(e(i)) for i in (1, 2, 3, 4)]
v1, u1 = bps_pair_random(); v2, u2 = bps_pair_random()
report('generic two-block BPS (companion)', [v1, u1, v2, u2])
report('generic global BPS', [v1, u1, v1, u1])
report('block1 generic, block2 vacuum', [v1, u1, e1, e4])
report('block1 vacuum, block2 generic', [e1, e4, v2, u2])
# flags: v in span(e1..e_{a}), u in span(e_{5-b}..e4) with v.u=0
def flag_pair(da, db):
    v = np.zeros(4, dtype=np.int64); u = np.zeros(4, dtype=np.int64)
    v[:da + 1] = rv(da + 1); u[3 - db:] = rv(db + 1)
    s = int((v * u % p).sum() % p)
    if s:
        # fix: adjust the u component with largest overlap index
        k = max(i for i in range(4) if v[i] and u[i])
        u[k] = (u[k] - s * pow(int(v[k]), -1, p)) % p
    return v, u
for da in range(4):
    for db in range(4):
        if da + db > 4: continue
        va, ua = flag_pair(da, db); vb, ub = flag_pair(da, db)
        try:
            report(f'flag deg(p)={da},deg(q)={db}: block1 only', [va, ua, e1, e4])
            report(f'flag deg(p)={da},deg(q)={db}: block2 only', [e1, e4, vb, ub])
            report(f'flag deg(p)={da},deg(q)={db}: both blocks (independent)', [va, ua, vb, ub])
        except Exception as ex:
            print(ex)
