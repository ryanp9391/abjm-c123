from common import *
# Q-tilde (twisted polynomial gauge) with monic q_j; the formula's Gamma factors are state independent
def Qt(s, j, x): return z[j]**(x / h) * qj(s, j, x)
def Q12(s, x): return Qt(s, 0, x) * Qt(s, 1, x - h) - Qt(s, 1, x) * Qt(s, 0, x - h)
def Q123(s, x): return det([[Qt(s, i, x - k * h) for k in range(3)] for i in range(3)])
for s in states[:12]:
    form = Qt(s, 0, th2) * Q12(s, th2) * Q123(s, th2) * Qt(s, 0, th1) * Q12(s, th1) * Q123(s, th1 - h)
    vac = s['ketm'][3]     # canonical ket p_a (x) V_b has <0|Psi> = 1
    # also the 'sigma-permuted' variants Q_1 -> Q_b etc. are not used
    print(s['label'], '<0|Psi>/formula =', mp.nstr(vac / form, 10))
