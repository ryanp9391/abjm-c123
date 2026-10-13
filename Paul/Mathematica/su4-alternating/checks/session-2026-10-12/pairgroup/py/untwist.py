# Untwisting test of the pair word at the ABJM point (python numerics, complex128)
#   twist z = 1 + eps*phi with phi2 = phi3 (z2 = z3), pair-preserving inhomogeneity theta_F,n = eps*thdir_n, theta_A = theta_F + 2h
#   mode 'vac'  : kappa = 1/z4, a = -1:  <vac (x) W_-1^{L-1}|Psi> = K(tau) <W_-1^L|Psi>/(1 - kappa z1)^2  -> <vac W_-1^{L-1}|Psi0>
#   mode 'shift': a = 1/3, kappa s.t. a' - a = -1 at every eps (kappa -> 1):  <W_a' (x) W_a^{L-1}|Psi> -> <vac^{L-1} W_-1|Psi0>
import numpy as np, itertools, sys
from pgbase import Chain, prod_state
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
N = 2*L; h = 0.5; D = 4**N
types = ['A', 'F']*L
allidx = np.array(list(itertools.product(range(4), repeat=N)))
wts = np.zeros((D, 4), dtype=int)
for s in range(N):
    sg = 1 if types[s] == 'A' else -1
    wts[np.arange(D), allidx[:, s]] += sg
sectors = {}
for n in range(D): sectors.setdefault(tuple(wts[n]), []).append(n)
perm_shift = np.array([np.ravel_multi_index(tuple(idx[2:]) + tuple(idx[:2]), (4,)*N) for idx in allidx])

def apply_shift(psi):
    return psi[perm_shift]

def apply_J(i, j, psi):
    """ket action of the global generator rho(E_ij) = sum_s rho_s(E_ij)"""
    T = psi.reshape((4,)*N); out = np.zeros_like(T)
    for s in range(N):
        m = np.zeros((4, 4))
        if types[s] == 'A': m[i, j] = 1
        else: m[j, i] = -1
        out += np.moveaxis(np.tensordot(m, T, axes=([1], [s])), 0, s)
    return out.reshape(-1)

def sector_matrix(ch, k, u0, sec, reg=False):
    """M[p,q] = (e_p . t_k(u0))_q restricted to the sector; ket action (t Psi) = M Psi."""
    d = len(sec); M = np.zeros((d, d), dtype=complex)
    for p, n in enumerate(sec):
        e = np.zeros(D, dtype=complex); e[n] = 1
        row = ch.apply_t(k, u0, e.reshape((4,)*N), reg=reg).reshape(-1)
        M[p] = row[sec]
    return M

def nullspace(A, tol=1e-10):
    u, s, vh = np.linalg.svd(A)
    r = np.sum(s > tol*max(1.0, s[0] if len(s) else 1.0))
    return vh[r:].conj().T

def chain(eps, phi, thdir):
    thF = [eps*thdir[n] for n in range(L)]
    th = []
    for f in thF: th += [f + 2*h, f]
    z = 1 + eps*np.array(phi)
    G = np.diag(z).astype(complex).astype(object)
    return Chain(types, th, h, G), th, z

ch0, th0, _ = chain(0.0, [0, 0, 0, 0], [0]*L)

def dominant(w):
    return all(w[i] >= w[i+1] for i in range(3))

hw_states = []
for w, sec in sectors.items():
    if not dominant(w): continue
    d = len(sec)
    rows = []
    for (i, j) in ((0, 1), (1, 2), (2, 3)):
        Jm = np.zeros((D, d))
        for p, n in enumerate(sec):
            e = np.zeros(D); e[n] = 1
            Jm[:, p] = apply_J(i, j, e)
        rows.append(Jm[np.abs(Jm).sum(1) > 0])
    Um = np.zeros((d, d))
    for p, n in enumerate(sec):
        e = np.zeros(D); e[n] = 1
        Um[:, p] = apply_shift(e)[sec]
    rows.append(Um - np.eye(d))
    NS = nullspace(np.vstack(rows))
    if NS.shape[1] == 0: continue
    Mg = sector_matrix(ch0, 1, 0.37, sec) + 0.61*sector_matrix(ch0, 3, -0.29, sec) + 0.23*sector_matrix(ch0, 2, 0.83, sec)
    R = np.linalg.lstsq(NS, Mg @ NS, rcond=None)[0]
    ev, V = np.linalg.eig(R)
    for m in range(len(ev)):
        psi = np.zeros(D, dtype=complex); psi[sec] = NS @ V[:, m]
        psi /= np.linalg.norm(psi)
        k = np.argmax(np.abs(psi)); psi *= np.abs(psi[k])/psi[k]
        hw_states.append((w, psi))
print('L =', L, ': number of HW cyclic joint eigenstates:', len(hw_states))
e = np.eye(4)

def W(a): return [e[0] + a*e[3], e[3] - a*e[0]]

def bra(vecs): return prod_state(vecs, dtype=complex).reshape(-1)

vac = [e[0], e[3]]
B_first = bra(vac + W(-1)*(L-1))      # vac (x) W_-1^{L-1}
B_last = bra(vac*(L-1) + W(-1))       # vac^{L-1} (x) W_-1
for w, psi in hw_states:
    a1 = B_first @ psi; a2 = B_last @ psi
    print('  weight', w, ' <vac W-1^{L-1}|Psi0> = %+.6f%+.6fi   <vac^{L-1} W-1|Psi0> = %+.6f%+.6fi' % (a1.real, a1.imag, a2.real, a2.imag))
np.save('hw_L%d.npy' % L, np.array([psi for w, psi in hw_states]))

def lam_hat(th, k, u0):
    """localisation constant of t_k(u0) on pair 1 (product over sites of pairs 2..L, exactly-vanishing factors dropped)"""
    r = 1.0
    for s in range(2, N):
        t = th[s]
        if s % 2 == 0: fac = [u0 - t - j*h for j in range(k)]
        else: fac = [u0 - t - k*h] + [u0 - t - j*h for j in range(k-1)]
        for f in fac:
            if f != 0: r *= f
    return r

def word_check(phi, thdir, eps, mode, a):
    ch, th, z = chain(eps, phi, thdir)
    y = z[1]
    u1 = th[1]; u2 = th[1] + h
    l1, l2, l3 = lam_hat(th, 1, u1), lam_hat(th, 2, u2), lam_hat(th, 3, u2)
    out = []
    for w, psi0 in hw_states:
        sec = sectors[w]
        Mg = sector_matrix(ch, 1, 0.37, sec) + 0.61*sector_matrix(ch, 3, -0.29, sec) + 0.23*sector_matrix(ch, 2, 0.83, sec)
        ev, V = np.linalg.eig(Mg)
        p0 = psi0[sec]
        best = None; used = np.zeros(len(ev), bool)
        for m in range(len(ev)):
            if used[m]: continue
            cl = np.where(np.abs(ev - ev[m]) < 1e-7*max(1, abs(ev[m])))[0]; used[cl] = True
            Vc = V[:, cl]
            c = np.linalg.lstsq(Vc, p0, rcond=None)[0]; proj = Vc @ c
            if best is None or np.linalg.norm(proj) > best[0]: best = (np.linalg.norm(proj), proj, len(cl))
        ov, ps, mult = best
        psi = np.zeros(D, dtype=complex); psi[sec] = ps/np.linalg.norm(ps)
        psi *= np.exp(-1j*np.angle(np.vdot(psi0, psi)))
        def tauk(k, u0, reg=False):
            M = sector_matrix(ch, k, u0, sec, reg=reg); v = M @ psi[sec]
            t = np.vdot(psi[sec], v)/np.vdot(psi[sec], psi[sec])
            return t, np.linalg.norm(v - t*psi[sec])/max(abs(t), 1e-300)
        t1, r1 = tauk(1, u1); t2, r2 = tauk(2, u2); t3, r3 = tauk(3, u2, reg=True)
        kap = 1/z[3] if mode == 'vac' else 1/(z[0] - a*(z[0] - z[3]))
        s = 1 - kap*(z[0] + z[3])
        K = s*s + kap*s*(t1/(h**2*l1) + t2/(2*h**4*y*l2) + t3/(6*h**5*y**2*l3)) + kap**2*t3*t1/(24*h**7*y**2*l1*l3)
        ap = a*(1 - kap*z[3])/(1 - kap*z[0])
        ovW = bra(W(a)*L) @ psi
        Nk = (1 - kap*z[0])**2
        rhs = K*ovW/Nk
        lhs = bra(W(ap) + W(a)*(L-1)) @ psi
        out.append(dict(w=w, ov=ov, mult=mult, eigres=max(r1, r2, r3), ovW=ovW, KN=K/Nk, rhs=rhs, lhs=lhs, ap=ap, kap=kap))
    return out

if __name__ == '__main__':
    phi = [0.9, 0.35, 0.35, -1.1]; thdir = [0.0, 0.37, -0.52][:L]
    for mode, a in (('vac', -1.0), ('shift', 1/3)):
        print('==== mode', mode, ' a =', a)
        for eps in (1e-2, 1e-3, 1e-4):
            rows = word_check(phi, thdir, eps, mode, a)
            for (w, psi0), r in zip(hw_states, rows):
                tgt = (B_first if mode == 'vac' else B_last) @ psi0
                print("  eps=%.0e w=%s ov=%.6f mult=%d eigres=%.0e |<W^L|Psi>|=%.2e |K/N|=%.2e RHS=%+.8f%+.8fi |RHS-LHS|=%.1e target=%+.8f |RHS-target|=%.1e a'=%.4f"
                      % (eps, r['w'], r['ov'], r['mult'], r['eigres'], abs(r['ovW']), abs(r['KN']), r['rhs'].real, r['rhs'].imag,
                         abs(r['rhs'] - r['lhs']), tgt.real, abs(r['rhs'] - tgt), r['ap'].real))
