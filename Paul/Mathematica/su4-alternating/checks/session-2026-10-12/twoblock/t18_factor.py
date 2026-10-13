from sympy import factorint, Rational
from fractions import Fraction as Fr
from core import *
import itertools
th = th0; h = h0; z = z0
c = esym(z)
print('e_k(z) =', c, [ (factorint(x.numerator), factorint(x.denominator)) for x in c])
print('theta diffs:')
for i, j in itertools.combinations(range(4), 2):
    for m in (-2, -1, 0, 1, 2):
        d = th[i] - th[j] + m * h
        print(f'  th{i+1}{j+1}{m:+d}h = {d}  ', factorint(d.numerator), factorint(d.denominator))
vals = [Fr(-2207744, 6933615), Fr(-174076198912, 739615645665), Fr(-14000, 2263), Fr(5390000, 1894131), Fr(-241472, 50895), Fr(-1822147712, 41805153), Fr(64,343), Fr(-10,21), Fr(-640,7203)]
for v in vals:
    print(v, factorint(v.numerator), factorint(v.denominator))
