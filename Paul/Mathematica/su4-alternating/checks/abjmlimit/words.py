import numpy as np
from fractions import Fraction as Fr
from itertools import product
from chain import *

H = Fr(1, 2)
Z = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
G = companion(Z)
TYPES = ['A', 'F', 'A', 'F']
WIDX = list(product(range(4), repeat=4))  # (a1,k2,a3,k4) lexicographic like Mathematica Tuples


def words(ch, uA1=None, uA3=None, ops=None):
    """256 words <0| t_a1(uA1) t1(th2)^k2 t_a3(uA3) t1(th4)^k4 (rows).
    ops: optional dict overriding the A-site operator lists: ops['A1'][a], ops['A3'][a]."""
    F = ch.F
    th = ch.th
    if uA1 is None:
        uA1 = th[0] - ch.h
    if uA3 is None:
        uA3 = th[2] - ch.h
    vac = vacuum(ch.types, F)
    A1 = ops['A1'] if ops and 'A1' in ops else [ch.t(a, uA1) for a in range(4)]
    A3 = ops['A3'] if ops and 'A3' in ops else [ch.t(a, uA3) for a in range(4)]
    T2 = ch.t(1, th[1])
    T4 = ch.t(1, th[3])
    rows = {}
    for a1 in range(4):
        v1 = F.red(vac @ A1[a1]) if isinstance(F, ModP) else vac @ A1[a1]
        v = v1
        for k2 in range(4):
            if k2 > 0:
                v = F.red(v @ T2) if isinstance(F, ModP) else v @ T2
            for a3 in range(4):
                w = F.red(v @ A3[a3]) if isinstance(F, ModP) else v @ A3[a3]
                for k4 in range(4):
                    if k4 > 0:
                        w = F.red(w @ T4) if isinstance(F, ModP) else w @ T4
                    rows[(a1, k2, a3, k4)] = w
    return np.array([rows[k] for k in WIDX])
