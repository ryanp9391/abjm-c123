"""Transfer matrices t1,t2,t3 of the alternating su(4) chain A F A F ... acting on bras (row vectors).

Conventions (task spec): L_A = (u-th) + h K, L_F = (u-th) - h P on aux x site, K = sum e_ij x e_ij, P = sum e_ij x e_ji.
Monodromy T = L_1 ... L_n (aux product), twist on the right.  Bras psi.M : new_b = sum_a psi_a M_ab.
Site 1 = most significant index of the 4^n vector.
dtype may be float64 or object (Fraction / mpmath mpf) for exact / high precision checks.
"""
import numpy as np
import itertools
from fractions import Fraction

NN = 4
H = Fraction(1, 2)
TH6 = [Fraction(2, 7), Fraction(-3, 5), Fraction(8, 11), Fraction(-11, 8), Fraction(5, 13), Fraction(-7, 9)]
TH6b = [Fraction(1, 3), Fraction(-2, 5), Fraction(3, 4), Fraction(-5, 4), Fraction(2, 9), Fraction(-7, 10)]
TH6c = [Fraction(19, 23), Fraction(34, 29), Fraction(41, 40), Fraction(-19, 12), Fraction(-7, 11), Fraction(13, 17)]
Z4 = [Fraction(3, 2), Fraction(5, 4), Fraction(2, 3), Fraction(4, 7)]
SUBS2 = list(itertools.combinations(range(4), 2))


def e_(i, j, dt=float):
    m = np.zeros((4, 4), dtype=dt)
    m[i, j] = 1
    return m


def rho_tensor(typ, dt):
    """R[k,l,a,b] = rho_typ(E_kl)[a,b]; A: e_kl, F: -e_lk."""
    R = np.zeros((4, 4, 4, 4), dtype=dt)
    for k in range(4):
        for l in range(4):
            if typ == 'A':
                R[k, l, k, l] = 1
            else:
                R[k, l, l, k] = -1
    return R


def wedgeD():
    """D(e_ab)_{I,J} as in dual-lax-t3.wb: Dw[X][Ir,Jc] = sum_k X[k,J1] wc[k,J2,Ir] + X[k,J2] wc[J1,k,Ir]."""
    def wc(k, m, Ir):
        if (k, m) == Ir:
            return 1
        if (m, k) == Ir:
            return -1
        return 0
    D = np.zeros((4, 4, 6, 6), dtype=int)
    for a in range(4):
        for b in range(4):
            X = np.zeros((4, 4), dtype=int)
            X[a, b] = 1
            for iI, Ir in enumerate(SUBS2):
                for iJ, Jc in enumerate(SUBS2):
                    D[a, b, iI, iJ] = sum(X[k, Jc[0]] * wc(k, Jc[1], Ir) + X[k, Jc[1]] * wc(Jc[0], k, Ir) for k in range(4))
    return D


DW = wedgeD()


def to_dt(x, dt):
    if dt is float:
        return float(x)
    return x


def lax_tensor(which, typ, u, th, h, dt=float):
    """Lax[k,l,a,b] (aux k,l; site a,b), as a dtype-dt array.  which in 1,2,3.  Also returns site prefactor g."""
    one = np.eye(4, dtype=dt)
    if which == 1:
        R = rho_tensor(typ, dt)
        d = 4
        Lx = np.zeros((4, 4, 4, 4), dtype=dt)
        for k in range(4):
            for l in range(4):
                if k == l:
                    Lx[k, l] = (u - th) * one
                Lx[k, l] = Lx[k, l] + h * R[k, l]
        g = 1
    elif which == 3:
        Lx = np.zeros((4, 4, 4, 4), dtype=dt)
        for k in range(4):
            for l in range(4):
                if typ == 'A':
                    if k == l:
                        Lx[k, l] = (u - th + h) * one
                    Lx[k, l] = Lx[k, l] - h * e_(l, k, dt)
                else:
                    if k == l:
                        Lx[k, l] = (u - th - 3 * h) * one
                    Lx[k, l] = Lx[k, l] + h * e_(k, l, dt)
        g = (u - th - h) * (u - th - 2 * h) if typ == 'A' else (u - th) * (u - th - h)
    elif which == 2:
        R = rho_tensor(typ, dt)
        Lx = np.zeros((6, 6, 4, 4), dtype=dt)
        shift = 0 if typ == 'A' else h
        for iI in range(6):
            for iJ in range(6):
                if iI == iJ:
                    Lx[iI, iJ] = (u - th - shift) * one
                acc = np.zeros((4, 4), dtype=dt)
                for a in range(4):
                    for b in range(4):
                        c = int(DW[a, b, iI, iJ])
                        if c != 0:
                            acc = acc + c * R[a, b]
                Lx[iI, iJ] = Lx[iI, iJ] + h * acc
        g = (u - th - h) if typ == 'A' else (u - th)
    return Lx, g


def twist_matrix(which, G):
    """twist matrix Gt with t = sum T_ij Gt_ji."""
    if which == 1:
        return G
    if which == 3:
        # det G * G^{-T}: cofactor-type, computed with numpy/object
        return cof(G)
    if which == 2:
        W = np.zeros((6, 6), dtype=G.dtype)
        for jJ, J in enumerate(SUBS2):
            for iI, I in enumerate(SUBS2):
                W[jJ, iI] = det_(G[np.ix_(J, I)])
        return W


def det_(M):
    n = M.shape[0]
    if n == 1:
        return M[0, 0]
    if n == 2:
        return M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    raise ValueError


def cof(G):
    """det G * G^{-T} = adj(G)^T, i.e. Gbar[l,k] = (-1)^{k+l} det G_{without row l? ...}. Use generic exact formula.
    (Gbar)_{lk} = (-1)^{k+l} det G_{cK, cL}... we compute adj(G)^T-like via minors:  Gbar = det(G) * inv(G)^T.
    """
    n = 4
    Gb = np.zeros((4, 4), dtype=G.dtype)
    for l in range(4):
        for k in range(4):
            rows = [r for r in range(4) if r != l]
            cols = [c for c in range(4) if c != k]
            sub = G[np.ix_(rows, cols)]
            Gb[l, k] = (-1) ** (k + l) * det3(sub)
    return Gb


def det3(M):
    return (M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
            - M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
            + M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]))


class Chain:
    def __init__(self, nsites, th, h, G, dt=float):
        self.n = nsites
        self.types = ['A' if i % 2 == 0 else 'F' for i in range(nsites)]
        self.th = th[:nsites]
        self.h = h
        self.dt = dt
        self.G = np.array(G, dtype=dt if dt is not float else float)
        self.tw = {w: twist_matrix(w, self.G) for w in (1, 2, 3)}
        self.dim = 4 ** nsites
        self._lax = {}

    def conv(self, x):
        return float(x) if self.dt is float else x

    def gprod(self, which, u):
        """prefactor prod_s g_s(u) computed in the native number type of th,h,u (exact if Fractions)."""
        g = 1
        for s in range(self.n):
            t = self.th[s]; h = self.h
            if which == 1:
                gs = 1
            elif which == 2:
                gs = (u - t - h) if self.types[s] == 'A' else (u - t)
            else:
                gs = (u - t - h) * (u - t - 2 * h) if self.types[s] == 'A' else (u - t) * (u - t - h)
            g = g * gs
        return self.conv(g)

    def laxes(self, which, u):
        key = (which, u)
        if key not in self._lax:
            if len(self._lax) > 400:
                self._lax.clear()
            Ls = []
            for s in range(self.n):
                Lx, g = lax_tensor(which, self.types[s], self.conv(u), self.conv(self.th[s]), self.conv(self.h), self.dt)
                Ls.append(Lx)
            gtot = self.gprod(which, u)
            self._lax[key] = (Ls, gtot)
        return self._lax[key]

    def apply(self, which, u, psi, prefactor=True):
        """psi: array (dim,) or (dim, B).  returns psi . t_which(u)."""
        Ls, gtot = self.laxes(which, u)
        d = 4 if which in (1, 3) else 6
        n = self.n
        batch = psi.ndim == 2
        B = psi.shape[1] if batch else 1
        ps = psi.reshape(self.dim, B)
        # Phi[i,k,state,B]
        Phi = np.zeros((d, d, self.dim, B), dtype=ps.dtype if self.dt is float else object)
        for i in range(d):
            Phi[i, i] = ps
        for s in range(n):
            Lx = Ls[s]
            Phi5 = Phi.reshape(d, d, 4 ** s, 4, 4 ** (n - s - 1) * B)
            # new[i,l,x,b,y] = sum_{k,a} Phi5[i,k,x,a,y] Lx[k,l,a,b]
            new = np.tensordot(Phi5, Lx, axes=([1, 3], [0, 2]))  # (i,x,y,l,b)
            new = np.transpose(new, (0, 3, 1, 4, 2))  # (i,l,x,b,y)
            Phi = np.ascontiguousarray(new).reshape(d, d, self.dim, B)
        Gt = self.tw[which]
        out = np.zeros((self.dim, B), dtype=Phi.dtype)
        for i in range(d):
            for j in range(d):
                c = Gt[j, i]
                if c != 0:
                    out = out + c * Phi[i, j]
        if prefactor:
            out = out * gtot if self.dt is float else out * gtot
        return out if batch else out[:, 0]


def build_chain(nsites, which_th=1, h=None, z=None, dt=float):
    th = TH6 if which_th == 1 else TH6b
    h = H if h is None else h
    z = Z4 if z is None else z
    if dt is float:
        G = np.diag([float(x) for x in z])
    else:
        G = np.diag(np.array(z, dtype=object)) if False else np.zeros((4, 4), dtype=object)
        for i in range(4):
            for j in range(4):
                G[i, j] = z[i] if i == j else (Fraction(0) if dt is Fraction else dt(0))
    return Chain(nsites, th, h, G, dt)


def site_vec(a):
    v = np.zeros(4)
    v[0] = 1.0
    v[3] = a
    u = np.zeros(4)
    u[3] = 1.0
    u[0] = -a
    return v, u


def kron_all(vs):
    out = vs[0]
    for v in vs[1:]:
        out = np.kron(out, v)
    return out


def W_state(nsites, a, a1=None):
    """W_{a1} on pair 1 (if given) otherwise W_a, times W_a^(pairs-1)."""
    vs = []
    for p in range(nsites // 2):
        aa = a1 if (p == 0 and a1 is not None) else a
        v, u = site_vec(aa)
        vs += [v, u]
    return kron_all(vs)
