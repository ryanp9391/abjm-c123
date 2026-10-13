"""Dual numbers a + b eps (eps^2 = 0) over Fractions, for exact u-derivatives of transfer matrices with pgbase."""
from fractions import Fraction as Fr
import numpy as np
class D:
    __slots__ = ('a', 'b')
    def __init__(self, a, b=0): self.a, self.b = a, b
    def _c(self, o): return o if isinstance(o, D) else D(o, 0)
    __array_priority__ = -1
    def __add__(self, o):
        if isinstance(o, np.ndarray): return NotImplemented
        o = self._c(o); return D(self.a + o.a, self.b + o.b)
    __radd__ = __add__
    def __sub__(self, o):
        if isinstance(o, np.ndarray): return NotImplemented
        o = self._c(o); return D(self.a - o.a, self.b - o.b)
    def __rsub__(self, o):
        if isinstance(o, np.ndarray): return NotImplemented
        o = self._c(o); return D(o.a - self.a, o.b - self.b)
    def __mul__(self, o):
        if isinstance(o, np.ndarray): return NotImplemented
        o = self._c(o); return D(self.a * o.a, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__
    def __neg__(self): return D(-self.a, -self.b)
    def __truediv__(self, o):
        if isinstance(o, np.ndarray): return NotImplemented
        o = self._c(o); return D(self.a / o.a, (self.b * o.a - self.a * o.b) / (o.a * o.a))
    def __rtruediv__(self, o):
        if isinstance(o, np.ndarray): return NotImplemented
        return D(o, 0) / self
    def __eq__(self, o): o = self._c(o); return self.a == o.a and self.b == o.b
    def __ne__(self, o): return not self.__eq__(o)
    def __hash__(self): return hash((self.a, self.b))
    def __repr__(self): return f'D({self.a},{self.b})'
def val(x): return x.a if isinstance(x, D) else x
def der(x): return x.b if isinstance(x, D) else 0
