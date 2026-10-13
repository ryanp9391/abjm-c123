import numpy as np, itertools
import eng
from eng import P
def support_pairs(ch, vec, alist):
    """smallest set S of pairs (0-based) with vec = (anything on S) (x) (W_{a_n} on pairs not in S)  (mod p)"""
    Lp = ch.Lp
    for r in range(0, Lp + 1):
        for S in itertools.combinations(range(Lp), r):
            V = vec.reshape((16,) * Lp) % P; ok = True
            for n in range(Lp):
                if n in S: continue
                w = eng.prod_state(eng.W(alist[n]), 'p').reshape(16)
                piv = int(np.nonzero(w)[0][0])
                Vn = np.moveaxis(V, n, -1)
                if not np.array_equal((Vn * w[piv]) % P, (Vn[..., piv:piv + 1] * w[None, :]) % P): ok = False; break
            if ok: return S
    return None
