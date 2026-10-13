import pickle
from taylor import *
d = pickle.load(open('regL4_both.pkl','rb')); R, Rd = d['R'], d['Rd']; p = P
def is_reg(a,k): return a == 3 or (a == 2 and k <= 2)
ratios = {}
for i,w in enumerate(WIDX):
    j = np.nonzero(R[i])[0][0]
    r = int(Rd[i][j]) * pow(int(R[i][j]), p-2, p) % p
    key = (w[0] if is_reg(w[0],w[1]) else '-', w[2] if is_reg(w[2],w[3]) else '-')
    ratios.setdefault(key, set()).add(ratrec(r))
for k,v in sorted(ratios.items(), key=str): print('reg a1,a3 =', k, 'ratio d-deriv/u-deriv:', v)
