"""ue_base.py -- epsilon-expansion machinery for the even-bridge identity (session 2026-10-13).
Reuses lpairs/eng.py (fused t1,t2,t3 on bras, diagonal twist) and lpairs/lp.py; works for Lp = 2, 3 pairs.
Conventions: chain A F A F ..., bras row vectors, psi -> psi.M; ket action of the same matrix = M psi.
Untwisted ABJM point: theta_F,n = eps*thdir_n, theta_A = theta_F + 2h, z = 1 + eps*phi (phi2 = phi3), h = 1/2.
HW cyclic states Psi0 from data/abjm-hw-states.wxf (keys L2, L3 only)."""
import numpy as np, itertools, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LP = os.path.normpath(os.path.join(HERE, '..', '..', 'session-2026-10-12', 'lpairs'))
sys.path.insert(0, LP)
import eng, lp, wxf
DATA = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'data', 'abjm-hw-states.wxf'))
h = 0.5
GEN = [(1.0, 1, 0.37, False), (0.61, 3, -0.29, True), (0.23, 2, 0.83, False)]   # generic commuting combination (as untw.py)

class Setup:
    def __init__(self, Lp):
        assert Lp in (2, 3)
        self.Lp = Lp; self.N = N = 2 * Lp; self.D = D = 4 ** N
        self.types = ['A', 'F'] * Lp
        tup = np.array(list(itertools.product(range(4), repeat=N)))
        wt = np.zeros((D, 4), dtype=int)
        for s in range(N):
            sg = 1 if self.types[s] == 'A' else -1
            for i in range(4): wt[:, i] += sg * (tup[:, s] == i)
        self.wt = wt
        sect = {}
        for n in range(D): sect.setdefault(tuple(wt[n]), []).append(n)
        self.SECT = {k: np.array(v) for k, v in sect.items()}
        self.load_states()
        # global generator rho(E14) as a sparse-ish matrix acting on kets (J psi); bra action psi.J is the same matrix
        self.J14 = self.global_gen(0, 3)
        self.J41 = self.global_gen(3, 0)

    def load_states(self):
        cache = os.path.join(HERE, 'hw_L%d.npz' % self.Lp)
        if os.path.exists(cache):
            d = np.load(cache, allow_pickle=True)
            self.PSI0, self.W0, self.DYN, self.E0, self.C0 = d['psi'], d['w'], d['dyn'], d['E'], d['C']
            return
        d = wxf.load(DATA)['L%d' % self.Lp]
        self.PSI0 = np.array([np.array(s['psi'], dtype=complex) for s in d])
        self.W0 = np.array([np.array(s['w'], dtype=int) for s in d])
        self.DYN = np.array([np.array(s['dyn'], dtype=int) for s in d])
        self.E0 = np.array([float(s['E']) for s in d])
        self.C0 = np.array([np.array(s['C'], dtype=complex) for s in d])
        np.savez(cache, psi=self.PSI0, w=self.W0, dyn=self.DYN, E=self.E0, C=self.C0)

    def global_gen(self, i, j):
        """ket-action matrix of rho(E_ij) = sum_s rho_s(E_ij) (A: e_ij, F: -e_ji) on the full space (dense for L=2, sparse for L=3)"""
        import scipy.sparse as sp
        N = self.N
        m_A = sp.csr_matrix(([1.0], ([i], [j])), shape=(4, 4))
        m_F = sp.csr_matrix(([-1.0], ([j], [i])), shape=(4, 4))
        I4 = sp.identity(4, format='csr')
        tot = None
        for s in range(N):
            fac = [I4] * N; fac[s] = m_A if self.types[s] == 'A' else m_F
            M = fac[0]
            for f in fac[1:]: M = sp.kron(M, f, format='csr')
            tot = M if tot is None else tot + M
        return tot.tocsr()

    def chain(self, eps, phi, thdir):
        thF = [eps * t for t in thdir[:self.Lp]]
        z = list(1 + eps * np.array(phi, dtype=float))
        return lp.chain_abjm(self.Lp, thF, h, z, 'c'), z

    def sector_matrix(self, ch, ops, sec, chunk=64):
        """M[p,q] = (e_p . sum_i c_i t_{k_i}(u_i))_q restricted to the sector; ket action M psi.  ops = [(c,k,u,reg)]"""
        d = len(sec); N = self.N; D = self.D; M = np.zeros((d, d), dtype=complex)
        for b0 in range(0, d, chunk):
            idx = sec[b0:b0 + chunk]; B = len(idx)
            E = np.zeros((D, B), dtype=complex); E[idx, np.arange(B)] = 1
            E = E.reshape((4,) * N + (B,))
            acc = np.zeros((D, B), dtype=complex)
            for c, k, u, reg in ops:
                acc += c * ch.apply_t(k, u, E, reg=reg).reshape(D, B)
            M[b0:b0 + B] = acc[sec].T
        return M

    def twisted_states(self, ch, which=None):
        """dict i -> (psi_eps full normalised phase-aligned, overlap, sector, gap)"""
        out = {}
        which = range(len(self.PSI0)) if which is None else which
        bysec = {}
        for i in which: bysec.setdefault(tuple(self.W0[i]), []).append(i)
        for w, ids in bysec.items():
            sec = self.SECT[w]
            Mg = self.sector_matrix(ch, GEN, sec)
            ev, V = np.linalg.eig(Mg)
            for i in ids:
                p0 = self.PSI0[i][sec]
                ovs = np.abs(V.conj().T @ p0) / np.linalg.norm(V, axis=0)
                m = int(np.argmax(ovs))
                v = V[:, m] / np.linalg.norm(V[:, m])
                v *= np.exp(-1j * np.angle(np.vdot(p0, v)))
                gaps = np.abs(ev - ev[m]); gaps[m] = np.inf
                psi = np.zeros(self.D, dtype=complex); psi[sec] = v
                out[i] = (psi, float(ovs[m] / np.linalg.norm(p0)), sec, float(gaps.min()))
        return out

    def eigval(self, ch, psi, sec, k, u, reg=False):
        M = self.sector_matrix(ch, [(1.0, k, u, reg)], sec)
        v = psi[sec]; Mv = M @ v
        t = np.vdot(v, Mv) / np.vdot(v, v)
        return t, float(np.linalg.norm(Mv - t * v) / max(abs(t) * np.linalg.norm(v), 1e-300))

    def pair_taus(self, ch, psi, sec, p=1):
        """(tau_F, tau2(u0), tau3^(u0)) = eigenvalues of X1, X2, X3 on psi, pair p; also x = tau1(u0)"""
        uF = ch.th[2 * p - 1]
        t1, r1 = self.eigval(ch, psi, sec, 1, uF)
        t2, r2 = self.eigval(ch, psi, sec, 2, uF + h)
        t3, r3 = self.eigval(ch, psi, sec, 3, uF + h, reg=True)
        x, r4 = self.eigval(ch, psi, sec, 1, uF + h)
        return (t1, t2, t3, x), max(r1, r2, r3, r4)

    @staticmethod
    def Kval(ch, taus, kap, y, p=1):
        c = lp.Kcoef(ch, p, kap, y)
        t1, t2, t3 = taus[:3]
        return c[0] + c[1] * t1 + c[2] * t2 + c[3] * t3 + c[4] * t3 * t1

    @staticmethod
    def Kcoef(ch, kap, y, p=1):
        return lp.Kcoef(ch, p, kap, y)

    def bra(self, alist):
        vecs = []
        for a in alist: vecs += eng.W(a)
        return eng.prod_state(vecs, 'c').reshape(-1)

    def vac(self):
        return eng.prod_state(eng.VAC * self.Lp, 'c').reshape(-1)

    def vacJn(self, n):
        """bra <vac| rho(E14)^n  (row vector)"""
        v = self.vac()
        for _ in range(n): v = self.J14.T @ v          # psi.J  == J^T psi as a column
        return v

def untwisted_chain(Lp):
    """exactly homogeneous untwisted ABJM chain (z = 1, theta_F = 0, theta_A = 2h); reduced t3 drops ALL vanishing factors"""
    return lp.chain_abjm(Lp, [0.0] * Lp, h, [1.0] * 4, 'c')

def lam_red(ch, k, u, p=1):
    return eng.lam(ch, k, u, p)

def K1_untwisted(ch, taus, p=1):
    """K at kappa = 1, z = 1 from the pair-point eigenvalues (tau_F, tau2(u0), tau3^(u0)) with reduced lambdas"""
    return Setup.Kval(ch, taus, 1.0, 1.0, p)

def nval(Lp, dyn):
    """number of E41 lowering steps from the vacuum HW weight to [p,0,p]: n = L - p (None if w2 or w3 != 0)"""
    p, q, r = dyn
    if q != 0 or p != r: return None
    return Lp - p
