import numpy as np, itertools
h = 0.5
def types(Lp): return ['A','F']*Lp
def _site_move(X, s, src, dst, coef, out):
    """out[..b=dst..] += coef * X[..a=src..]  along tensor axis s (axis index s in X, where X may have leading aux axis -> caller passes axis)"""
    sl_src=[slice(None)]*X.ndim; sl_src[s]=src
    sl_dst=[slice(None)]*X.ndim; sl_dst[s]=dst
    out[tuple(sl_dst)] += coef * X[tuple(sl_src)]
def t1_apply(psi, th, z, u, h=h, bar=False, N=None):
    """psi: array shape (4,)*N (+ optional batch axes). returns psi . t(u). bar=False: t1 ; bar=True: tbar3 (dual Lax, twist prod(z)/z)."""
    N = len(th)
    ty = types(N//2)
    tw = np.array(z) if not bar else np.prod(z)/np.array(z)
    res = np.zeros(psi.shape, dtype=np.result_type(psi.dtype, np.array(z).dtype, float))
    for i in range(4):
        X = np.zeros((4,)+psi.shape, dtype=res.dtype); X[i]=psi
        for s in range(N):
            x = u - th[s]
            Xn = np.zeros_like(X)
            for k in range(4):
                for kp in range(4):
                    # site op M_{k,kp}
                    if k==kp:
                        c = x + (h if (bar and ty[s]=='A') else 0) + (-3*h if (bar and ty[s]=='F') else 0)
                        Xn[kp] += c * X[k]
                    if not bar:
                        if ty[s]=='A': _site_move(X[k], s, k, kp, h, Xn[kp])      # h e_{k kp}: a=k -> b=kp
                        else: _site_move(X[k], s, kp, k, -h, Xn[kp])             # -h e_{kp k}: a=kp -> b=k
                    else:
                        if ty[s]=='A': _site_move(X[k], s, kp, k, -h, Xn[kp])     # -h e_{kp k}
                        else: _site_move(X[k], s, k, kp, h, Xn[kp])               # +h e_{k kp}
            X = Xn
        res += tw[i]*X[i]
    return res
def E(i): return np.eye(4)[i-1]
def Wv(a): return [E(1)+a*E(4), E(4)-a*E(1)]
def prod_state(sv):
    r=np.array([1.0])
    for v in sv: r=np.kron(r,v)
    return r.reshape((4,)*len(sv))
def rot(sv, s, M):
    sv=[np.array(v) for v in sv]; sv[s]=sv[s]@M; return sv
def ratio(x, y, tol=1e-13):
    xf=x.reshape(-1); yf=y.reshape(-1)
    k=np.abs(yf)>tol*np.abs(yf).max(); r=xf[k]/yf[k]
    return r[np.argmax(np.abs(yf[k]))], np.max(np.abs(r-r[np.argmax(np.abs(yf[k]))]))/abs(r[0]), (np.max(np.abs(xf[~k])) if (~k).any() else 0.0)
def thetas(Lp, th): return list(th[:2*Lp])
TH8 = [2/7, -3/5, 8/11, -11/8, 5/13, -7/9, 4/15, -9/11]
Z4 = [3/2, 5/4, 2/3, 4/7]
def vartheta(ty, k, s, th, forward):  # factor
    pass
# closed forms --------------------------------------------------
def d_closed(k, th, Lp, h=h):
    """step constant d_k (k=1-based): odd (A): <P_k| t1(th_k) = (1/d_k) <P_{k-1}| ; even (F): <P_{k-1}| t1(th_k) = d_k <P_k|"""
    N=2*Lp; ty=types(Lp)
    if k%2==1:
        p=h
        for s in range(1,N+1):
            if s==k: continue
            p*= (th[k-1]-th[s-1]) if ty[s-1]=='F' else (th[k-1]-th[s-1]+h)
        return 1/p
    else:
        p=-h
        for s in range(1,N+1):
            if s==k: continue
            p*= (th[k-1]-th[s-1]) if ty[s-1]=='A' else (th[k-1]-th[s-1]-h)
        return p
def c_closed(k, th, Lp):
    r=1.0
    for j in range(1,k+1): r*=d_closed(j,th,Lp)
    return r
