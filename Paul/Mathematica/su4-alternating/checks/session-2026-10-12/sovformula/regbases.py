"""Regularised SoV bases at the ABJM point (pair by pair, thF1 != thF2), L = 2 pairs (4 sites), exact Fractions.
Bra (ABJM-point agent): <0| O1 t1(th2)^k2 O3 t1(th4)^k4,  O in {1, t1(u0), W = t2'(u0) - t1'(u0) t1(thF) (k<=2) | t2(u0) (k=3), Z = t3'(u0)}.
Ket (mirror, this file): t1(th1)^m1 K2 t1(th3)^m3 K4 |0bar>,  K in {1, t1(u0), Wt = t2'(thA) - t1(thA) t1'(u0), Zt = t3'(thA+h)}
(at the ABJM point t2(thA) = t1(thA) t1(u0) and t3(thA+h) = 0 exactly, mirror of the bra identities)."""
import numpy as np, itertools, pickle, os
from fractions import Fraction as Fr
import sov_l4 as S
from dual import D, val, der

def dsplit(M):
    V = np.vectorize(val, otypes=[object])(M); Dm = np.vectorize(der, otypes=[object])(M)
    return V, Dm

class Reg:
    def __init__(self, h=Fr(1, 2), thF=(Fr(-3, 5), Fr(-11, 8)), z=S.Z0):
        self.h = h; self.thF = thF
        self.th = (thF[0] + 2 * h, thF[0], thF[1] + 2 * h, thF[1])
        self.z = z
        self.ch = S.Chain(h=h, th=self.th, z=z)
        self.chd = self.ch          # only the spectral parameter is dual
        self._ops = {}
    def op(self, a, u):
        return self.ch.t(a, u)
    def dop(self, a, u):
        """(t_a(u), t_a'(u)) exact"""
        key = (a, u)
        if key not in self._ops:
            M = self.chd.t(a, D(u, Fr(1)))
            self._ops[key] = dsplit(M)
        return self._ops[key]
    def pair_ops_bra(self, p):
        th = self.th; uF = th[2 * p + 1]; u0 = uF + self.h
        t1F = self.op(1, uF); t1u0, d1 = self.dop(1, u0); t2u0, d2 = self.dop(2, u0); t3u0, d3 = self.dop(3, u0)
        assert all(x == 0 for x in t3u0.reshape(-1))
        W = d2 - d1.dot(t1F)
        return dict(t1F=t1F, O=[None, t1u0, W, d3], t2u0=t2u0)
    def pair_ops_ket(self, p):
        th = self.th; thA = th[2 * p]; uF = th[2 * p + 1]; u0 = uF + self.h
        t1A, d1A = self.dop(1, thA); t1u0 = self.op(1, u0); t2A, d2A = self.dop(2, thA); t3, d3 = self.dop(3, thA + self.h)
        assert all(x == 0 for x in t3.reshape(-1))
        t1u0v, d1u0 = self.dop(1, u0)
        Wt = d2A - t1A.dot(d1u0)          # rank-256 choice (t14_ketvariants.py: 'Wb')
        return dict(t1A=t1A, K=[None, t1u0, Wt, d3])
    def bra_words(self):
        vac = S.sov_vacuum()
        P = [self.pair_ops_bra(0), self.pair_ops_bra(1)]
        rows = []
        for (a1, k2, a3, k4) in S.WORDS:
            v = vac
            for p, (a, k) in enumerate(((a1, k2), (a3, k4))):
                if a:
                    Op = P[p]['O'][a] if not (a == 2 and k == 3) else P[p]['t2u0']
                    v = v.dot(Op)
                for _ in range(k): v = v.dot(P[p]['t1F'])
            rows.append(v)
        return np.array(rows)
    def ket_words(self):
        vac = S.product_vec([S.unit(i) for i in S.KET_VAC])
        P = [self.pair_ops_ket(0), self.pair_ops_ket(1)]
        cols = []
        for (m1, n2, m3, n4) in S.WORDS:
            v = vac
            for p, (m, n) in reversed(list(enumerate(((m1, n2), (m3, n4))))):
                if n: v = P[p]['K'][n].dot(v)
                for _ in range(m): v = P[p]['t1A'].dot(v)
            cols.append(v)
        return np.array(cols).T
