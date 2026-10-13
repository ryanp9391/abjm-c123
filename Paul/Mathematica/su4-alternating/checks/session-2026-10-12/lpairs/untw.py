"""untw.py -- twisted eigenstates near the untwisted homogeneous ABJM chain (L = 3 pairs) and the pair-word formula.
Twisted chain: z = 1 + eps*phi, theta_F,n = eps*thdir_n, theta_A,n = theta_F,n + 2h (pair-preserving), h = 1/2.
Untwisted HW cyclic states Psi0 from data/abjm-hw-states.wxf (L3).  Twisted state = eigenvector (generic combination of
t's in the weight sector of Psi0) with maximal overlap with Psi0; normalised, phase aligned."""
import numpy as np, itertools
import eng, lp, hwdata
Lp = 3; N = 2 * Lp; D = 4 ** N; h = 0.5
types = ['A', 'F'] * Lp
_tup = np.array(list(itertools.product(range(4), repeat=N)))
_wt = np.zeros((D, 4), dtype=int)
for s in range(N):
    sg = 1 if types[s] == 'A' else -1
    for i in range(4): _wt[:, i] += sg * (_tup[:, s] == i)
SECT = {}
for n in range(D): SECT.setdefault(tuple(_wt[n]), []).append(n)
SECT = {k: np.array(v) for k, v in SECT.items()}
PSI0, W0, DYN, E0, C0 = hwdata.load()

def make_chain(eps, phi, thdir, zfun=None):
    thF = [eps * t for t in thdir]
    z = list(1 + eps * np.array(phi, dtype=float)) if zfun is None else zfun(eps)
    return lp.chain_abjm(Lp, thF, h, z, 'c'), z

def sector_matrix(ch, ops, sec, chunk=48):
    """sum_i c_i t_{k_i}(u_i) restricted to the sector, ket action M psi (M[p,q] = (e_p . t)_q). ops = [(c,k,u,reg)]"""
    d = len(sec); M = np.zeros((d, d), dtype=complex)
    for b0 in range(0, d, chunk):
        idx = sec[b0:b0 + chunk]; B = len(idx)
        E = np.zeros((D, B), dtype=complex); E[idx, np.arange(B)] = 1
        E = E.reshape((4,) * N + (B,))
        acc = np.zeros((D, B), dtype=complex)
        for c, k, u, reg in ops:
            acc += c * ch.apply_t(k, u, E, reg=reg).reshape(D, B)
        M[b0:b0 + B] = acc[sec].T
    return M

GEN = [(1.0, 1, 0.37, False), (0.61, 3, -0.29, True), (0.23, 2, 0.83, False)]

def twisted_states(ch, which=None):
    """returns dict state-index -> (psi_eps (full, normalised, phase-aligned), overlap, sector, eigen-gap)"""
    out = {}
    which = range(len(PSI0)) if which is None else which
    bysec = {}
    for i in which: bysec.setdefault(tuple(W0[i]), []).append(i)
    for w, ids in bysec.items():
        sec = SECT[w]
        Mg = sector_matrix(ch, GEN, sec)
        ev, V = np.linalg.eig(Mg)
        for i in ids:
            p0 = PSI0[i][sec]
            ovs = np.abs(V.conj().T @ p0) / np.linalg.norm(V, axis=0)
            m = int(np.argmax(ovs))
            v = V[:, m] / np.linalg.norm(V[:, m])
            v *= np.exp(-1j * np.angle(np.vdot(p0, v)))
            gaps = np.abs(ev - ev[m]); gaps[m] = np.inf
            psi = np.zeros(D, dtype=complex); psi[sec] = v
            out[i] = (psi, float(ovs[m] / np.linalg.norm(p0)), sec, float(gaps.min()))
    return out

def eigval(ch, psi, sec, k, u, reg=False):
    """eigenvalue of t_k(u) on the ket psi (+ relative residual)"""
    M = sector_matrix(ch, [(1.0, k, u, reg)], sec)
    v = psi[sec]; Mv = M @ v
    t = np.vdot(v, Mv) / np.vdot(v, v)
    return t, float(np.linalg.norm(Mv - t * v) / max(abs(t) * np.linalg.norm(v), 1e-300))

def pair_taus(ch, psi, sec, p=1):
    uF = ch.th[2 * p - 1]
    t1, r1 = eigval(ch, psi, sec, 1, uF)
    t2, r2 = eigval(ch, psi, sec, 2, uF + h)
    t3, r3 = eigval(ch, psi, sec, 3, uF + h, reg=True)
    return (t1, t2, t3), max(r1, r2, r3)

def Kval(ch, taus, kap, y, p=1):
    c = lp.Kcoef(ch, p, kap, y)
    t1, t2, t3 = taus
    return c[0] + c[1] * t1 + c[2] * t2 + c[3] * t3 + c[4] * t3 * t1

def bra(alist):
    vecs = []
    for a in alist: vecs += eng.W(a)
    return eng.prod_state(vecs, 'c').reshape(-1)

def weight_factor(z, kap, w):
    """global torus rho(1 - kappa Lambda) eigenvalue on a ket of weight w"""
    return np.prod([(1 - kap * z[i]) ** w[i] for i in range(4)])
