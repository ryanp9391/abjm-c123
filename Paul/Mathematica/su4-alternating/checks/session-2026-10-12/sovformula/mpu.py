"""High-precision (mpmath) untwisting: exact rational chain (Fractions) at rational eps, sector matrices of transfer
matrices (bra action restricted to a weight sector, ket action T[p,q]), joint eigenvectors with mpmath.
Use: sec = Sector(L, eps, phi, thdir); sec.build(w, ops); psi, psiL = sec.eig(w, psi0_sector)."""
import numpy as np, itertools, mpmath as mp
from fractions import Fraction as Fr
from pgbase import Chain, prod_state, diagm
H = Fr(1, 2)

def tomp(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x)

class Ctx:
    def __init__(self, L, eps, phi, thdir, dps=60):
        mp.mp.dps = dps
        self.L, self.N = L, 2 * L
        self.types = ['A', 'F'] * L
        self.eps = Fr(eps); self.phi = [Fr(p) for p in phi]
        thF = [self.eps * Fr(t) for t in thdir]
        self.th = []
        for f in thF: self.th += [f + 2 * H, f]
        self.z = [1 + self.eps * p for p in self.phi]
        self.ch = Chain(self.types, self.th, H, diagm(self.z))
        allidx = list(itertools.product(range(4), repeat=self.N))
        self.allidx = allidx
        self.sectors = {}
        for n, idx in enumerate(allidx):
            w = [0] * 4
            for s, i in enumerate(idx): w[i] += 1 if self.types[s] == 'A' else -1
            self.sectors.setdefault(tuple(w), []).append(n)
        z = self.z
        self.lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]

    def secmat(self, w, k, u, reg=False):
        """exact Fraction matrix M[p,q] = <e_p| t_k(u) |e_q>, p,q in sector w"""
        sec = self.sectors[w]; pos = {n: i for i, n in enumerate(sec)}
        M = [[Fr(0)] * len(sec) for _ in sec]
        for i, n in enumerate(sec):
            e = np.empty((4,) * self.N, dtype=object); e[...] = Fr(0)
            e[self.allidx[n]] = Fr(1)
            r = self.ch.apply_t(k, u, e, reg=reg).reshape(-1)
            for m in sec:
                if r[m] != 0: M[i][pos[m]] = r[m]
        return M

    def eig(self, w, Mgen, psi0):
        """joint eigenvector of the sector matrix Mgen (ket action) closest to psi0 (numpy complex, full space);
        returns (psi, psiL) mp column/row, <psiL|psi> = 1, psi normalised & phase aligned with psi0, overlap, gap"""
        sec = self.sectors[w]
        A = mp.matrix([[tomp(x) for x in row] for row in Mgen])
        E, ER = mp.eig(A)
        p0 = [mp.mpc(complex(psi0[n])) for n in sec]
        best, bm = -1, None
        for m in range(len(E)):
            v = [ER[i, m] for i in range(len(sec))]
            nv = mp.sqrt(mp.fsum(abs(x) ** 2 for x in v))
            o = abs(mp.fsum(mp.conj(a) * b for a, b in zip(p0, v))) / nv
            if o > best: best, bm = o, m
        gap = min([abs(E[m] - E[bm]) for m in range(len(E)) if m != bm] + [mp.inf])
        EL = mp.inverse(ER)
        psi = [ER[i, bm] for i in range(len(sec))]
        psiL = [EL[bm, i] for i in range(len(sec))]
        nv = mp.sqrt(mp.fsum(abs(x) ** 2 for x in psi))
        ph = mp.fsum(mp.conj(a) * b for a, b in zip(p0, psi)); ph = ph / abs(ph)
        psi = [x / (nv * ph) for x in psi]; psiL = [x * nv * ph for x in psiL]
        return psi, psiL, best, gap

    def prodbra(self, vecs, w):
        """components of the product bra (list of per-site vectors, Fractions/mp) on sector w"""
        out = []
        for n in self.sectors[w]:
            idx = self.allidx[n]; p = mp.mpf(1)
            for s, i in enumerate(idx):
                x = vecs[s][i]
                p *= tomp(x) if isinstance(x, (Fr, int)) else x
            out.append(p)
        return out

def dot(a, b): return mp.fsum(x * y for x, y in zip(a, b))
def matvec(M, v): return [mp.fsum(tomp(M[i][j]) * v[j] for j in range(len(v)) if M[i][j] != 0) for i in range(len(M))]
def vecmat(v, M): return [mp.fsum(v[i] * tomp(M[i][j]) for i in range(len(v)) if M[i][j] != 0) for j in range(len(M[0]))]
