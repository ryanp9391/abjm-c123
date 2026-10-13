import numpy as np, itertools, time
from scipy.optimize import linprog


def omp(A, t, kmax, tol=1e-9):
    """orthogonal matching pursuit on unit-norm columns of A (r x N). Returns selected indices, coefficients, residual history."""
    nr = np.linalg.norm(A, axis=0)
    nr[nr == 0] = 1
    An = A / nr
    res = t.copy()
    sel = []
    hist = []
    tn = np.linalg.norm(t)
    for k in range(kmax):
        corr = np.abs(An.T @ res)
        corr[sel] = 0
        j = int(np.argmax(corr))
        sel.append(j)
        c, *_ = np.linalg.lstsq(An[:, sel], t, rcond=None)
        res = t - An[:, sel] @ c
        hist.append(np.linalg.norm(res) / tn)
        if hist[-1] < tol:
            break
    c = c / nr[sel]
    return sel, c, hist


def l1(A, t, scale_cols=True, tol=1e-9):
    """min sum |x_j| subject to A x = t via HiGHS (x = p - n). Columns scaled to unit norm.
    Returns support, coefficients (in original column scaling), success flag."""
    r, N = A.shape
    nr = np.linalg.norm(A, axis=0) if scale_cols else np.ones(N)
    nr[nr == 0] = 1
    An = A / nr
    Aeq = np.concatenate([An, -An], axis=1)
    cost = np.ones(2 * N)
    res = linprog(cost, A_eq=Aeq, b_eq=t, bounds=(0, None), method='highs')
    if res.status != 0:
        return None, None, False
    x = res.x[:N] - res.x[N:]
    x[np.abs(x) < 1e-10 * np.abs(x).max()] = 0
    sup = np.nonzero(x)[0]
    return sup, x[sup] / nr[sup], True


def subset_search(A, t, k, cand=None, tol=1e-9, maxhits=20, tmax=600):
    """exhaustive k-subset search: subsets of columns of A (r x N) whose span contains t (rel. residual < tol)."""
    r, N = A.shape
    cand = list(range(N)) if cand is None else list(cand)
    nr = np.linalg.norm(A, axis=0)
    nr[nr == 0] = 1
    An = A / nr
    tn = np.linalg.norm(t)
    hits = []
    t0 = time.time()
    for sub in itertools.combinations(cand, k):
        M = An[:, sub]
        c, *_ = np.linalg.lstsq(M, t, rcond=None)
        rr = np.linalg.norm(M @ c - t) / tn
        if rr < tol:
            # require all coefficients nonzero (minimal)
            if np.all(np.abs(c) > 1e-9 * np.abs(c).max()):
                hits.append((sub, c / nr[list(sub)], rr))
                if len(hits) >= maxhits:
                    break
        if time.time() - t0 > tmax:
            hits.append(('timeout',))
            break
    return hits
