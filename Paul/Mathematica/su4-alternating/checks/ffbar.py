import numpy as np, itertools
N=4; rng=np.random.default_rng(1)
def e(i,j,n=N):
    m=np.zeros((n,n)); m[i,j]=1; return m
rho_f  = lambda i,j: e(i,j)
rho_fb = lambda i,j: -e(j,i)
h=0.5; th=[0.2857,-0.6]
def lax(rho,u):  # aux x phys : sum e_ij (x) (u d_ij + h rho(E_ij)), returned as 4x4 block of phys matrices
    return [[ (u*(i==j))*np.eye(N)+h*rho(i,j) for j in range(N)] for i in range(N)]
reps=[rho_f,rho_fb]
I=np.eye(N)
def T(u):  # T_ij = sum_k L1_ik L2_kj, ops on site1 (x) site2
    L1=lax(reps[0],u-th[0]); L2=lax(reps[1],u-th[1])
    return [[sum(np.kron(L1[i][k],L2[k][j]) for k in range(N)) for j in range(N)] for i in range(N)]
z=np.array([1.5,1.25,2/3,4/7])
chi=[np.poly(-z)[k] for k in range(N+1)]  # prod(t+z_i) coefficients: chi_0=1,...
G=np.zeros((N,N)); 
for j in range(N): G[0,j]=(-1)**j*chi[j+1]
for i in range(1,N): G[i,i-1]=1
print("eig G",np.sort(np.linalg.eigvals(G).real))
def t1(u):
    Tm=T(u); return sum(Tm[i][j]*G[j,i] for i in range(N) for j in range(N))
def rhoG(rho,g):  # group action: for f, g ; for fbar, (g^-1)^T
    return g if rho is rho_f else np.linalg.inv(g).T
K=sum(np.kron(e(i,j),e(i,j)) for i in range(N) for j in range(N))
P=sum(np.kron(e(i,j),e(j,i)) for i in range(N) for j in range(N))
for a in (0,1):
    M=t1(th[a]); print("site",a+1,"rank t1(theta)",np.linalg.matrix_rank(M,1e-9))
    # try M = (c0 + c1 X) (rho(G) on site a) with X in {K,P}, either order
    for name,R in (("K",K),("P",P)):
      for gname,g in (("G",G),("G^T",G.T),("G^-1",np.linalg.inv(G)),("G^-T",np.linalg.inv(G).T)):
        Ga=np.kron(g,I) if a==0 else np.kron(I,g)
        for order in ("RG","GR"):
            basis=[Ga, (R@Ga if order=="RG" else Ga@R)]
            A=np.array([b.ravel() for b in basis]).T
            c,res,*_=np.linalg.lstsq(A,M.ravel(),rcond=None)
            r=np.linalg.norm(A@c-M.ravel())/np.linalg.norm(M)
            if r<1e-9: print("  t1(theta_%d) = %.6g*%s_%d + %.6g*%s.%s  [%s]"%(a+1,c[0],gname,a+1,c[1],name,gname,order))
