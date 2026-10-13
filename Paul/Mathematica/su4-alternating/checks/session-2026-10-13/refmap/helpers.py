from rm_base import *
def prod_(xs):
    p = Fr(1)
    for x in xs: p *= x
    return p
def tw_u(w, z):
    return [prod_([z[j] - w[m] for m in range(3)]) for j in range(4)]
def tw_form(w, z):
    v = []
    for j in range(4):
        p = Fr(1)
        for m in range(4):
            if m != j: p *= (w[0] - z[m]) / (z[j] - z[m])
        v.append(p)
    return v, tw_u(w, z)
def comp_v(z):
    out = []
    for j in range(4):
        p = Fr(1)
        for m in range(4):
            if m != j: p *= (z[j] - z[m])
        out.append(1 / p)
    return out
def pred(v, u):
    """closed-form admissibility criterion at z2 = z3: blockwise BPS (v1u1 = v4u4 = 0, v2u2+v3u3 = 0) or support in {1,4}"""
    return (v[0] * u[0] == 0 and v[3] * u[3] == 0 and v[1] * u[1] + v[2] * u[2] == 0) or (v[1] == v[2] == u[1] == u[2] == 0)
