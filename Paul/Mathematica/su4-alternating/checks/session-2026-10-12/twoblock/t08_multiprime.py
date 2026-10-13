import numpy as np, pickle, time
from sympy import prevprime
from core import *
t0 = time.time()
primes = []
q = 2 ** 26
while len(primes) < 8:
    q = prevprime(q); primes.append(q)
res = {}
for p in primes:
    R = Ring('mod', p); ch = Chain(R)
    M, vac = sov_words(ch)
    Y = R.inv(M)
    res[p] = (M, Y)
    print(p, time.time() - t0, flush=True)
pickle.dump(res, open('dual_mod8.pkl', 'wb'))
