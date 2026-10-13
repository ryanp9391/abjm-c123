from fractions import Fraction as Fr
def nullspace(rows, n):
    """right nullspace of a list of rows (Fractions), n columns. returns list of basis vectors"""
    A = [list(r) for r in rows if any(x != 0 for x in r)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == len(A): break
    free = [c for c in range(n) if c not in piv]
    out = []
    for fc in free:
        v = [Fr(0)] * n; v[fc] = Fr(1)
        for i, pc in enumerate(piv): v[pc] = -A[i][fc]
        out.append(v)
    return out
def rank(rows, n): return n - len(nullspace(rows, n)) if rows else 0
