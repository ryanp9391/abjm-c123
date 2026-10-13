"""Untwisting machinery (numpy complex128), L pairs A F A F ..., h = 1/2, ABJM point pair by pair:
thF_n = eps*thdir_n, thA_n = thF_n + 2h, diagonal twist z = 1 + eps*phi.  Untwisted HW cyclic eigenstates as in pairgroup/untwist.py.
Provides: dense operators t_k(u) (bra action psi -> psi @ T), weight sectors, right/left twisted eigenvectors matched to Psi0,
companion <-> diagonal frame maps, SoV vacua <0| = (lambda (x) 1)^L, |0bar> = (1 (x) lambda)^L (diag frame)."""
import numpy as np, itertools
from pgbase import Chain, prod_state
h = 0.5

class Setup:
    def __init__(self, L):
        self.L = L; self.N = N = 2 * L; self.D = D = 4 ** N
        self.types = ['A', 'F'] * L
        allidx = np.array(list(itertools.product(range(4), repeat=N)))
        wts = np.zeros((D, 4), dtype=int)
        for s in range(N):
            sg = 1 if self.types[s] == 'A' else -1
            wts[np.arange(D), allidx[:, s]] += sg
        self.wts = wts
        self.sectors = {}
        for n in range(D): self.sectors.setdefault(tuple(wts[n]), []).append(n)
        self.perm_shift = np.array([np.ravel_multi_index(tuple(idx[2:]) + tuple(idx[:2]), (4,) * N) for idx in allidx])
        self.hw = self._hw_states()

    def chain(self, eps, phi, thdir):
        thF = [eps * thdir[n] for n in range(self.L)]
        th = []
        for f in thF: th += [f + 2 * h, f]
        z = 1 + eps * np.array(phi, dtype=float)
        G = np.diag(z).astype(complex).astype(object)
        return Chain(self.types, th, h, G), th, z

    def op(self, ch, k, u, reg=False):
        """dense bra-action matrix T: (psi t_k(u)) = psi @ T"""
        D, N = self.D, self.N
        T = np.zeros((D, D), dtype=complex)
        for n in range(D):
            e = np.zeros(D, dtype=complex); e[n] = 1
            T[n] = ch.apply_t(k, u, e.reshape((4,) * N), reg=reg).reshape(-1)
        return T

    def opder(self, ch, k, u, reg=False, dl=0.05, npts=None):
        """dense d/du t_k(u) by exact polynomial interpolation (degree <= 4L+? in u)"""
        npts = npts or (4 * self.L + 3)
        us = [u + dl * (j - npts // 2) for j in range(npts)]
        Ts = [self.op(ch, k, x, reg=reg) for x in us]
        # Lagrange derivative weights at u
        w = []
        for j, xj in enumerate(us):
            # d/du of l_j(u)
            s = 0.0
            for m, xm in enumerate(us):
                if m == j: continue
                p = 1.0 / (xj - xm)
                for q, xq in enumerate(us):
                    if q in (j, m): continue
                    p *= (u - xq) / (xj - xq)
                s += p
            w.append(s)
        return sum(wj * Tj for wj, Tj in zip(w, Ts))

    def apply_J(self, i, j, psi):
        N = self.N
        T = psi.reshape((4,) * N); out = np.zeros_like(T)
        for s in range(N):
            m = np.zeros((4, 4))
            if self.types[s] == 'A': m[i, j] = 1
            else: m[j, i] = -1
            out += np.moveaxis(np.tensordot(m, T, axes=([1], [s])), 0, s)
        return out.reshape(-1)

    def _hw_states(self):
        D = self.D
        ch0, _, _ = self.chain(0.0, [0, 0, 0, 0], [0] * self.L)
        Mg = self.op(ch0, 1, 0.37) + 0.61 * self.op(ch0, 3, -0.29) + 0.23 * self.op(ch0, 2, 0.83)
        out = []
        for w, sec in self.sectors.items():
            if not all(w[i] >= w[i + 1] for i in range(3)): continue
            d = len(sec); rows = []
            for (i, j) in ((0, 1), (1, 2), (2, 3)):
                Jm = np.zeros((D, d))
                for p, n in enumerate(sec):
                    e = np.zeros(D); e[n] = 1
                    Jm[:, p] = self.apply_J(i, j, e)
                rows.append(Jm[np.abs(Jm).sum(1) > 0])
            Um = np.zeros((d, d))
            for p, n in enumerate(sec):
                e = np.zeros(D); e[n] = 1
                Um[:, p] = e[self.perm_shift][sec]
            rows.append(Um - np.eye(d))
            A = np.vstack(rows)
            u_, s_, vh = np.linalg.svd(A)
            r = np.sum(s_ > 1e-10 * max(1.0, s_[0]))
            NS = vh[r:].conj().T
            if NS.shape[1] == 0: continue
            Mk = Mg[np.ix_(sec, sec)]            # ket action on the sector: (T psi)
            R = np.linalg.lstsq(NS, Mk @ NS, rcond=None)[0]
            ev, V = np.linalg.eig(R)
            for m in range(len(ev)):
                psi = np.zeros(D, dtype=complex); psi[sec] = NS @ V[:, m]
                psi /= np.linalg.norm(psi)
                k = np.argmax(np.abs(psi)); psi *= np.abs(psi[k]) / psi[k]
                out.append((w, psi))
        return out

    def twisted_pair(self, Mbra, w, psi0):
        """right eigenket Psi (ket action Mbra^T) in sector w closest to psi0 (normalised, phase aligned with psi0) and
        left eigenbra PsiL (row vector, PsiL Mbra = mu PsiL), normalised <PsiL|Psi> = 1."""
        sec = self.sectors[w]
        Mk = Mbra[np.ix_(sec, sec)]
        ev, V = np.linalg.eig(Mk)
        p0 = psi0[sec]
        ov = np.abs(V.conj().T @ p0) / np.linalg.norm(V, axis=0)
        m = int(np.argmax(ov))
        gap = np.min(np.abs(np.delete(ev, m) - ev[m])) if len(ev) > 1 else np.inf
        Vi = np.linalg.inv(V)
        psi = np.zeros(self.D, dtype=complex); psi[sec] = V[:, m]
        psiL = np.zeros(self.D, dtype=complex); psiL[sec] = Vi[m]
        nrm = np.linalg.norm(psi); psi /= nrm; psiL *= nrm
        ph = np.vdot(psi0, psi); psi *= abs(ph) / ph; psiL *= ph / abs(ph)
        return psi, psiL, ov[m], gap

# frames
def lam_of(z):
    z = np.asarray(z)
    return np.array([1 / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)])

def vand(z):
    z = np.asarray(z); return np.array([[zj ** (3 - i) for zj in z] for i in range(4)])

def comp_to_diag_site(z):
    """A: v_d = v_c V^{-T};  F: u_d = u_c V"""
    V = vand(z); return np.linalg.inv(V).T, V

def nprime(z):
    lam = lam_of(z); r = lam / lam[0]
    n = np.array([[1, 1, 1, r[3]], [0, 1, 0, -r[1]], [0, 0, 1, -r[2]], [0, 0, 0, 1]], dtype=complex)
    Dp = np.diag([1, r[1], r[2], 1]).astype(complex)
    return n, Dp

def load_data_states(L):
    """saved untwisted HW cyclic states (data/abjm-hw-states.wxf): list of dict(w, dyn, E, C, psi)"""
    import wxf
    Dd = wxf.load(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf')
    out = []
    for s in Dd['L%d' % L]:
        sa = s['psi']
        dims = int(sa.args[1][0]); rowptr, cols = sa.args[3][1]; vals = sa.args[3][2]
        v = np.zeros(dims, dtype=complex)
        for c, x in zip(np.array(cols).ravel(), vals): v[int(c) - 1] = x
        out.append(dict(w=tuple(int(x) for x in s['w']), dyn=list(s['dyn']), E=s['E'], C=s['C'], psi=v))
    return out
