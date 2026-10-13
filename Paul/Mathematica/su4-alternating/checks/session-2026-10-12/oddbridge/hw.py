import numpy as np, itertools, wxf, time
from ob import *
DATA="C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf"
_d=None
def hwdata():
    global _d
    if _d is None: _d=wxf.load(DATA)
    return _d
def sectors(Lp):
    N=2*Lp; ty=types(Lp)
    tup=np.array(list(itertools.product(range(4),repeat=N)))
    w=np.zeros((len(tup),4),dtype=int)
    for s in range(N):
        sg=1 if ty[s]=='A' else -1
        for i in range(4): w[:,i]+=sg*(tup[:,s]==i)
    d={}
    for idx,row in enumerate(map(tuple,w)): d.setdefault(row,[]).append(idx)
    return {k:np.array(v) for k,v in d.items()}
def mixvec(psi, th, z):
    return (t1_apply(psi,th,z,0.3711)+0.37*t1_apply(psi,th,z,-0.917)+0.29*t1_apply(psi,th,z,0.41,bar=True)-0.13*t1_apply(psi,th,z,-0.77,bar=True))
def sector_matrix(Lp, idx, th, z):
    N=2*Lp; n=len(idx)
    B=np.zeros((4**N,n),dtype=complex); B[idx,np.arange(n)]=1
    out=mixvec(B.reshape((4,)*N+(n,)), th, z).reshape(4**N,n)
    return out[idx,:].T   # M_mat[i,j]=(e_i.M)_j ; ket action psi -> M_mat @ psi
