import numpy as np, itertools
from states import load
S=load()
def digits(L):
    n=2*L; N=4**n
    idx=np.arange(N)
    D=np.empty((N,n),dtype=np.int8)
    for s in range(n):
        D[:,s]=(idx//4**(n-1-s))%4
    return D
DIG={L:digits(L) for L in (2,3,4)}
def weight_of_basis(L):
    D=DIG[L]; W=np.zeros((D.shape[0],4),dtype=int)
    for s in range(2*L):
        sign=1 if s%2==0 else -1   # site 1 (index0) = A
        for i in range(4):
            W[:,i]+=sign*(D[:,s]==i)
    return W
def apply_E(psi,L,i,j):
    """global rho(E_ij) on ket psi: A: e_ij, F: -e_ji (matrix units acting on column kets)"""
    n=2*L; T=psi.reshape([4]*n); out=np.zeros_like(T)
    for s in range(n):
        # A: E_ij|e_j> = |e_i>; F: -e_ji : |e_i> -> -|e_j>
        if s%2==0: src,dst,sg=j,i,1
        else: src,dst,sg=i,j,-1
        Ts=np.moveaxis(T,s,0)
        Os=np.moveaxis(out,s,0)
        Os[dst]+=sg*Ts[src]
    return out.reshape(-1)
def overlap(psi,L,vecs):
    """vecs: list of 2L row vectors (bra site vectors)"""
    T=psi.reshape([4]*(2*L))
    for v in vecs[::-1]:
        T=T@np.asarray(v,dtype=complex)   # contracts last axis
    return complex(T)
if __name__=='__main__':
    for L in (2,3,4):
        W=weight_of_basis(L)
        worst=0; wworst=0
        for s in S[L]:
            psi=s['psi']
            nz=np.abs(psi)>1e-13
            ws=np.unique(W[nz],axis=0)
            assert len(ws)==1 and (ws[0]==s['w']).all(), (ws,s['w'])
            for i in range(4):
                for j in range(i+1,4):
                    worst=max(worst,np.abs(apply_E(psi,L,i,j)).max())
        print(L,'weights ok; max |E_ij psi| (i<j)=',worst)
    # C check
    e=np.eye(4)
    for L in (2,3):
        bad=0
        for s in S[L]:
            for l in range(2*L+1):
                vecs=[]
                for site in range(2*L):
                    A=(site%2==0)
                    if site<l: vecs.append(e[0] if A else e[3])
                    else: vecs.append(e[0]-e[3] if A else e[3]+e[0])
                c=overlap(s['psi'],L,vecs)
                bad=max(bad,abs(c-s['C'][l]))
        print(L,'C_l reproduced, max diff',bad)
