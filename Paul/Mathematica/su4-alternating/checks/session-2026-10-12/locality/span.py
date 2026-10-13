import numpy as np

def orth_basis(W, tol=1e-9, chunk=800, verbose=False):
    """Incremental rank-revealing orthonormal basis of the column span of W (columns normalised first).
    Returns Q (dim x r) and the list of residual singular values seen (for the gap diagnostics)."""
    dim, N = W.shape
    nr = np.linalg.norm(W, axis=0)
    keep = nr > 0
    Q = np.zeros((dim, 0))
    gaps = []
    for c0 in range(0, N, chunk):
        C = W[:, c0:c0 + chunk]
        n = np.linalg.norm(C, axis=0)
        C = C[:, n > 0] / n[n > 0]
        for _ in range(2):
            if Q.shape[1]:
                C = C - Q @ (Q.T @ C)
        U, s, _ = np.linalg.svd(C, full_matrices=False)
        gaps.append(s)
        r = int((s > tol).sum())
        if r:
            Qn = U[:, :r]
            # re-orthogonalise against Q
            Qn = Qn - Q @ (Q.T @ Qn)
            Qn, _ = np.linalg.qr(Qn)
            Q = np.concatenate([Q, Qn], axis=1)
        if verbose: print('chunk', c0, 'rank so far', Q.shape[1], flush=True)
    return Q, gaps

def resid(Q, t):
    t = np.asarray(t, float)
    r = t - Q @ (Q.T @ t)
    r = r - Q @ (Q.T @ r)
    return np.linalg.norm(r) / np.linalg.norm(t)

def localised(Q, w0, nloc):
    """intersection of span(Q) with C^16 (x) w0 : returns cos-angle spectrum, x-basis (16-vectors) for those ~1."""
    dim = Q.shape[0]
    R = dim // 16
    wh = w0 / np.linalg.norm(w0)
    Bl = np.kron(np.eye(16), wh[:, None])  # (16*R, 16)
    A = Q.T @ Bl
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    return s, U, Vt

def targets(L, a=1/3, z=(1.5, 1.25, 2/3, 4/7)):
    z = np.array(z)
    v = np.array([1, 0, 0, a]); u = np.array([-a, 0, 0, 1.0])
    rest = W_state_rest(L, a)
    T = {}
    for a1 in (-2/3, 0.0, 2/3, 5/4):
        vv = np.array([1, 0, 0, a1]); uu = np.array([-a1, 0, 0, 1.0])
        T[f'T1 a\'={a1:+.4g}'] = np.kron(np.kron(vv, uu), rest)
    for k in (1/3, -2/5, 5/3):
        D = 1 - k * z
        T[f'T2 kappa={k:+.4g}'] = np.kron(np.kron(v * D, u / D), rest)
    for k in (1/3, -2/5, 5/3):
        T[f'T3a v(1-k/L)xu k={k:+.4g}'] = np.kron(np.kron(v * (1 - k / z), u), rest)
    T['T3b vL^-1 x u(1+L)'] = np.kron(np.kron(v / z, u * (1 + z)), rest)
    T['T3b vL^-1 x u(2-L)'] = np.kron(np.kron(v / z, u * (2 - z)), rest)
    return T

from lib import site_vec, kron_all
def W_state_rest(L, a):
    vs = []
    for p in range(1, L // 2):
        v, u = site_vec(a); vs += [v, u]
    return kron_all(vs) if vs else np.array([1.0])
