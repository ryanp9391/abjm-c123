import numpy as np, functools
N=4; h=0.5; th=[0.2857,-0.6]
def e(i,j): m=np.zeros((N,N)); m[i,j]=1; return m
rf=lambda i,j: e(i,j); rb=lambda i,j: -e(j,i); I=np.eye(N)
def comp(z):
    c=np.poly(-np.array(z)); G=np.zeros((N,N)); G[0,:]=[(-1)**j*c[j+1] for j in range(N)]
    for i in range(1,N): G[i,i-1]=1
    return G
def Tm(u):
    L1=[[np.kron((u-th[0])*(i==j)*I+h*rf(i,j),I) for j in range(N)] for i in range(N)]
    L2=[[np.kron(I,(u-th[1])*(i==j)*I+h*rb(i,j)) for j in range(N)] for i in range(N)]
    return [[sum(L1[i][k]@L2[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
def ta(a,u,G):  # sum over I,J of qminor_{I,J}[T](u) * det G[J,I]
    Ts={}
    def T(v):
        if v not in Ts: Ts[v]=Tm(v)
        return Ts[v]
    out=np.zeros((16,16))
    for I_ in itertools.combinations(range(N),a):
        for J in itertools.combinations(range(N),a):
            q=np.zeros((16,16))
            for p in itertools.permutations(range(a)):
                sg=np.linalg.det(np.eye(a)[list(p)])
                q+=sg*functools.reduce(np.matmul,[T(u-k*h)[I_[p[k]]][J[k]] for k in range(a)])
            out+=q*np.linalg.det(G[np.ix_(J,I_)])
    return out
import itertools
vac=np.kron(np.eye(N)[0],np.eye(N)[3])
def basis(z):
    G=comp(z); B=[]
    t2=ta(1,th[1],G)
    for a in range(4):
        v=vac if a==0 else vac@ta(a,th[0]-h,G)
        for k in range(4):
            B.append(v@np.linalg.matrix_power(t2,k))
    return np.array(B)
B1=basis([1.5,1.25,2/3,4/7]); B2=basis([0.9,-1.7,2.2,3.1])
print("rank",np.linalg.matrix_rank(B1,1e-8)," max|B(z)-B(z')|/max|B| = %.1e"%(np.abs(B1-B2).max()/np.abs(B1).max()))
