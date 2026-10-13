import pickle
import mpmath as mp
from frames import *
we = pickle.load(open('words_exact.pkl', 'rb')); W = we['words']
ge = pickle.load(open('gram_exact.pkl', 'rb')); Y = ge['Y']
Beig = dict(pickle.load(open('Beig.pkl', 'rb')))
for n in Beig:
    k = companion_ket(n); b = companion_bra(n)
    ov0 = mp.fsum(tomp(W[0][c]) * k[c] for c in range(256) if W[0][c] != 0)
    ovL0 = mp.fsum(b[c] * tomp(Y[c, 0]) for c in range(256) if Y[c, 0] != 0)
    norm = mp.fsum(b[c] * k[c] for c in range(256))
    print(n, 'direct <PsiL|Psi>/(<0|Psi><PsiL|0bar>) =', mp.nstr(norm / (ov0 * ovL0), 15), ' from Q-measure:', mp.nstr(Beig[n], 15))
