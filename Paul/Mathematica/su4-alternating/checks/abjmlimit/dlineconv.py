import pickle, numpy as np
from fractions import Fraction as Fr
from words import WIDX
R = pickle.load(open('exact_R0.pkl','rb'))
def is_reg(a,k): return a == 3 or (a == 2 and k <= 2)
def pair_terms(a,k): return [((a,k),1)] if not (a==2 and k<=2) else [((2,k),1),((1,k+1),-1)]
for tag in ('100','1000','10000'):
    M = pickle.load(open('exact_M_%s.pkl'%tag,'rb')); delta = Fr(1,int(tag))
    worst = 0
    for j,w in enumerate(WIDX):
        n = is_reg(w[0],w[1]) + is_reg(w[2],w[3])
        row = np.zeros(256, dtype=object); row[:] = Fr(0)
        for (q1,c1) in pair_terms(w[0],w[1]):
            for (q3,c3) in pair_terms(w[2],w[3]):
                row = row + c1*c3*M[WIDX.index((q1[0],q1[1],q3[0],q3[1]))]
        row = row / delta**n
        dev = max(abs(float(x - y)) for x,y in zip(row, R[j])) / max(abs(float(y)) for y in R[j])
        worst = max(worst, dev)
    print('delta = 1/%s: max_j |(S M(delta))_j / delta^n_j - R_j| / |R_j| =' % tag, '%.3e' % worst)
