"""Task 5: word support of global Borel rotations of the companion SoV vacuum <0| = e1 e4 e1 e4 by g = 1 + b E_ij (diagonal-frame generator
transported to the companion frame: A: g_c = V^{-T} g V^T, F: V g^{-T} V^{-1})."""
from wa_sov_core import *
import numpy as np
ch, W, Y, Yk, G, M = build(verbose=False)
Yn = np.array(Y, dtype=object)
V = vandermonde(Z0); Vi = minv(V)
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
def tr(A): return [list(r) for r in zip(*A)]
def E(i, j): return [[Fr(1) if (r == i and c == j) else Fr(0) for c in range(4)] for r in range(4)]
Vt = tr(V); Vit = tr(Vi)
def rowvec_times(v, X): return [sum(v[k] * X[k][j] for k in range(4)) for j in range(4)]
for (i, j) in [(0, 1), (2, 3), (0, 3), (0, 2), (1, 3), (1, 2)]:
    XA = mm(mm(Vit, E(i, j)), Vt)                 # V^{-T} E V^T
    XF = [[-x for x in r] for r in mm(mm(V, E(j, i)), Vi)]     # -V E_ji V^{-1}
    e1, e4 = unit(0), unit(3)
    A0, A1 = e1, rowvec_times(e1, XA)
    F0, F1 = e4, rowvec_times(e4, XF)
    comp = {n: np.array([Fr(0)] * 256, dtype=object) for n in range(5)}
    for flips in itertools.product([0, 1], repeat=4):
        vs = [(A1 if fl else A0) if s % 2 == 0 else (F1 if fl else F0) for s, fl in enumerate(flips)]
        comp[sum(flips)] = comp[sum(flips)] + product_vec(vs)
    cs = {n: comp[n].dot(Yn) for n in comp}
    cnt = [sum(1 for x in cs[n] if x != 0) for n in range(5)]
    tot = sum(1 for k in range(256) if any(cs[n][k] != 0 for n in range(5)))
    print('E%d%d' % (i + 1, j + 1), 'words per b-power', cnt, 'union', tot, flush=True)

print('--- contrast: generators acting directly in the companion frame (g_c = 1 + b E_ij, A: v g_c, F: u g_c^{-T})')
for (i, j) in [(0, 1), (2, 3), (0, 3), (0, 2), (1, 3), (1, 2)]:
    XA = E(i, j); XF = [[-x for x in r] for r in E(j, i)]
    e1, e4 = unit(0), unit(3)
    A0, A1 = e1, rowvec_times(e1, XA)
    F0, F1 = e4, rowvec_times(e4, XF)
    comp = {n: np.array([Fr(0)] * 256, dtype=object) for n in range(5)}
    for flips in itertools.product([0, 1], repeat=4):
        vs = [(A1 if fl else A0) if s % 2 == 0 else (F1 if fl else F0) for s, fl in enumerate(flips)]
        comp[sum(flips)] = comp[sum(flips)] + product_vec(vs)
    cs = {n: comp[n].dot(Yn) for n in comp}
    cnt = [sum(1 for x in cs[n] if x != 0) for n in range(5)]
    tot = sum(1 for k in range(256) if any(cs[n][k] != 0 for n in range(5)))
    print('companion-frame E%d%d' % (i + 1, j + 1), 'words per b-power', cnt, 'union', tot, flush=True)
