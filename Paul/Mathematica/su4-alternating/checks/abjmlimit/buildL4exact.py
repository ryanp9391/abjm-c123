import pickle, time, sys
from taylor import *
from crt import exact_from_primes
h = H
tF1, tF2 = Fr(-3,5), Fr(-11,8)
USTAR = Fr(3137,10000)
def th_of(e): return [tF1+2*h+e, tF1, tF2+2*h+e, tF2]
def words_p(delta):
    def fn(p):
        ch = Chain(ModP(p), TYPES, th_of(delta), h, G); return words(ch)
    return fn
def t1_p(delta):
    def fn(p):
        ch = Chain(ModP(p), TYPES, th_of(delta), h, G); return ch.t(1, USTAR)
    return fn
def reg_p(p):
    F = ModP(p); th = th_of(Fr(0))
    ch = Chain(F, TYPES, th, h, G)
    TA1 = [t_taylor(ch, a, th[0]-h, 1) for a in range(4)]
    TA3 = [t_taylor(ch, a, th[2]-h, 1) for a in range(4)]
    Y2 = ch.t(1, th[1]); Y4 = ch.t(1, th[3])
    W1 = (TA1[2][1] - F.mm(TA1[1][1], Y2)) % p; Z1 = TA1[3][1]
    W3 = (TA3[2][1] - F.mm(TA3[1][1], Y4)) % p; Z3 = TA3[3][1]
    def op(TA, W, Z, a, k):
        if a < 2: return TA[a][0]
        if a == 2: return W if k <= 2 else TA[2][0]
        return Z
    vac = vacuum(TYPES, F); rows = []
    for (a1,k2,a3,k4) in WIDX:
        v = vac @ op(TA1,W1,Z1,a1,k2) % p
        for _ in range(k2): v = v @ Y2 % p
        v = v @ op(TA3,W3,Z3,a3,k4) % p
        for _ in range(k4): v = v @ Y4 % p
        rows.append(v)
    return np.array(rows)
job = sys.argv[1]
t0 = time.time()
if job == 'reg':
    R = exact_from_primes(reg_p, verbose=True); pickle.dump(R, open('exact_R0.pkl','wb'))
    T = exact_from_primes(t1_p(Fr(0)), verbose=True); pickle.dump(T, open('exact_t1_0.pkl','wb'))
else:
    delta = Fr(1, int(job))
    M = exact_from_primes(words_p(delta), verbose=True); pickle.dump(M, open('exact_M_%s.pkl'%job,'wb'))
    T = exact_from_primes(t1_p(delta), verbose=True); pickle.dump(T, open('exact_t1_%s.pkl'%job,'wb'))
print(job, 'done', time.time()-t0, flush=True)
