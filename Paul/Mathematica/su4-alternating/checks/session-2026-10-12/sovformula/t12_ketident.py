"""ABJM point (pair by pair), L = 2, exact: operator identities relevant for the ket SoV basis on |0bar> = e4 e1 e4 e1."""
import numpy as np
from fractions import Fraction as Fr
import sov_l4 as S
from dual import D, val, der
h = Fr(1, 2); thF = [Fr(-3, 5), Fr(-11, 8)]
th = (thF[0] + 2 * h, thF[0], thF[1] + 2 * h, thF[1])
ch = S.Chain(h=h, th=th, z=S.Z0)
thA1, thF1 = th[0], th[1]
u0 = thF1 + h
def eqz(M): return all(x == 0 for x in M.reshape(-1))
t1A = ch.t(1, thA1); t1u0 = ch.t(1, u0); t1F = ch.t(1, thF1)
print('t2(u0) == t1(u0) t1(thF)      :', eqz(ch.t(2, u0) - t1u0.dot(t1F)), ' (bra identity check, operators commute)')
print('t3(u0) == 0                  :', eqz(ch.t(3, u0)))
print('t2(thA) == t1(thA) t1(u0)    :', eqz(ch.t(2, thA1) - t1A.dot(t1u0)))
print('t3(thA+h) == 0               :', eqz(ch.t(3, thA1 + h)))
print('t2(thF+2h)=t2(thA): rank of t2(thA) - c t1(thA)t1(u0)? (diff zero?)', eqz(ch.t(2, thF1 + 2 * h) - ch.t(2, thA1)))
for (k, u, name) in [(2, thA1, 't2(thA)'), (3, thA1 + h, 't3(thA+h)'), (1, thA1, 't1(thA)')]:
    pass
