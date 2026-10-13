"""W_w(Psi) = <0_w|Psi_z>/(<0_w|0bar><0|Psi_z>) = sum_{x,y} tauTilde_y(vac_w) M_yx tau_x(Psi) on all 256 eigenstates (mp, 50 digits);
also the biorthogonality W_z(Psi) = 0 for Psi != vac, and sector structure of <W_a|Psi>."""
from wa_sov_core import *
import numpy as np
import sov_l4 as S
import mpmath as mp
mp.mp.dps = 50
ch, W, Y, Yk, G, M = build(verbose=False)
eig = cached('eigen', S.eigenstates)
sts = [S.Eig(d) for d in eig]
right, left = cached('lr', lambda: S.companion_vectors(eig))
print('eigenstates', len(sts), flush=True)
Ykn = np.array(Yk, dtype=object)
Wt = (Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9))
def vac_bra(zz):
    Vz = vandermonde(zz); Viz = minv(Vz)
    v = [Vz[i][0] for i in range(4)]; uu = Viz[3]
    return product_vec([v, uu, v, uu])
b0w = vac_bra(Wt); b0z = vac_bra(Z0)
vacsov = sov_vacuum()
idx = {w: i for i, w in enumerate(WORDS)}
def tomp_(x): return S.tomp(x)
def Wdirect(bra, n):
    num = mp.fsum(tomp_(bra[c]) * right[n][c] for c in range(256) if bra[c] != 0)
    den = mp.fsum(right[n][c] for c in range(256) if vacsov[c] != 0)   # <0| = e1 e4 e1 e4 : single component
    return num, den
b0wy = b0w.dot(Ykn); b0zy = b0z.dot(Ykn)
res = {}
errs = []; biorth = []
for n, st in enumerate(sts):
    tau_site = [[(mp.mpf(1) if lab == 0 else (st.tau(lab, S._th[s] - S._h) if S.TYPES[s] == 'A' else st.tau(1, S._th[s]) ** lab)) for lab in range(4)] for s in range(4)]
    tx = {x: mp.fprod([tau_site[s][x[s]] for s in range(4)]) for x in WORDS}
    for name, by in (('w', b0wy), ('z', b0zy)):
        tot = mp.fsum(tomp_(by[idx[y]] / by[0]) * tomp_(v) * tx[x] for (y, x), v in M.items())
        res[(name, n)] = tot
    num, den = Wdirect(b0w, n)
    direct = num / (tomp_(b0wy[0]) * den)
    errs.append(abs(res[('w', n)] / direct - 1))
    numz, _ = Wdirect(b0z, n)
    directz = numz / (tomp_(b0zy[0]) * den)
    biorth.append((st.weight, abs(res[('z', n)]), abs(directz)))
print('max_n |formula/direct - 1| for W_w:', max(errs))
vacn = [n for n, st in enumerate(sts) if st.weight == (2, 0, 0, -2)][0]
print('vac index', vacn, ' W_z(vac) formula:', res[('z', vacn)])
print('max_{n != vac} |W_z(Psi_n)| (formula):', max(b[1] for n, b in enumerate(biorth) if n != vacn), ' direct:', max(b[2] for n, b in enumerate(biorth) if n != vacn))
print('min |W_w(Psi_n)| over n:', min(abs(res[('w', n)]) for n in range(256)))
