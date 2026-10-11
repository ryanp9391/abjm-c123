import numpy as np, functools
N=4; h=0.5
def e(i,j): m=np.zeros((N,N)); m[i,j]=1; return m
rf=lambda i,j: e(i,j); rb=lambda i,j: -e(j,i)
z=np.array([1.5,1.25,2/3,4/7]); chi=np.poly(-z)
G=np.zeros((N,N)); G[0,:]=[(-1)**j*chi[j+1] for j in range(N)]
for i in range(1,N): G[i,i-1]=1
def chain(reps,th):
    L=len(reps); I=np.eye(N)
    def site_op(a,M): return functools.reduce(np.kron,[M if b==a else I for b in range(L)])
    def t1(u):
        # T = L1 L2 ... LL (aux product), entries ops
        T=[[np.eye(N**L)*(i==j) for j in range(N)] for i in range(N)]
        for a,r in enumerate(reps):
            La=[[site_op(a,(u-th[a])*(i==j)*I+h*r(i,j)) for j in range(N)] for i in range(N)]
            T=[[sum(T[i][k]@La[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
        return sum(T[i][j]*G[j,i] for i in range(N) for j in range(N))
    return t1,site_op
# L=2: diagonal frame check
t1,so=chain([rf,rb],[0.2857,-0.6])
M=t1(0.13); w,V=np.linalg.eig(M)
Pd=np.linalg.eig(G)  # G = W diag W^-1
print("L=2 t1 eigenvalues distinct:",len(np.unique(np.round(w,8))))
# rotate to diagonal frame: twist G=W D W^-1; act on kets: site1 by ? test product structure: count eigvecs of rank-1 as 4x4 matrix
ranks=[np.linalg.matrix_rank(V[:,k].reshape(N,N),1e-8) for k in range(16)]
print("L=2 eigvec matrix ranks:",sorted(ranks))
# L=4 alternating
th=[0.2857,-0.6,0.731,-1.37]
t1,so=chain([rf,rb,rf,rb],th)
I=np.eye(N)
K=lambda a,b: sum(so(a,e(i,j))@so(b,e(i,j)) for i in range(N) for j in range(N))
P=lambda a,b: sum(so(a,e(i,j))@so(b,e(j,i)) for i in range(N) for j in range(N))
for a in range(4):
    M=t1(th[a]); print("L=4 site",a+1,"rank",np.linalg.matrix_rank(M,1e-8))
# commutation: does site-pair (1,2) twist rotation commute with t1(u)?
R12=so(0,G.T)@so(1,np.linalg.inv(G))
print("[(G^T on sites1,2), t1(u)] norm:",np.linalg.norm(R12@t1(0.37)-t1(0.37)@R12))
M1=t1(th[0]); M2=t1(th[1])
X=M1@np.linalg.inv(M2)
print("|| t(th1) t(th2)^-1 - block rot ||/||X||:",np.linalg.norm(X-R12)/np.linalg.norm(X))
