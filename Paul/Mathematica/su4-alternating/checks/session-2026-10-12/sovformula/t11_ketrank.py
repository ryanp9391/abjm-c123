"""Rank of the ket SoV basis |m1,n2,m3,n4> = t1(th1)^m1 t_n2(th2+n2 h) t1(th3)^m3 t_n4(th4+n4 h)|0bar> at the ABJM point
(pair by pair, thF1 != thF2) and of the standard bra word basis there (exact, mod p)."""
import sys, numpy as np
from fractions import Fraction as Fr
import sov_l4 as S
h = Fr(1, 2); thF = [Fr(-3, 5), Fr(-11, 8)]
th = (thF[0] + 2 * h, thF[0], thF[1] + 2 * h, thF[1])
ch = S.Chain(h=h, th=th, z=S.Z0)
Y = S.ket_words(ch)
print('ABJM point: rank ket words (mod p) =', S.rank_mod_p(Y.T))
W = S.sov_words(ch)
print('ABJM point: rank standard bra words (mod p) =', S.rank_mod_p(W))
