import numpy as np, functools, itertools
N=4; h=0.5
def e(i,j): m=np.zeros((N,N)); m[i,j]=1; return m
rA=lambda i,j: e(i,j)      # notebook 'f' = Paul's anti-fund (built with t_n)
rF=lambda i,j: -e(j,i)     # notebook 'fbar' = Paul's fund (built with t1^n)
I=np.eye(N)
def comp(z):
    c=np.poly(-np.array(z)); G=np.zeros((N,N)); G[0,:]=[(-1)**j*c[j+1] for j in range(N)]
    for i in range(1,N): G[i,i-1]=1
    return G
G=comp([1.5,1.25,2/3,4/7])
class Chain:
    def __init__(s,reps,th): s.reps,s.th,s.L,s.cache=reps,th,len(reps),{}
    def so(s,a,M): return functools.reduce(np.kron,[M if b==a else I for b in range(s.L)])
    def T(s,u):
        if u in s.cache: return s.cache[u]
        D=N**s.L; T=[[np.eye(D)*(i==j) for j in range(N)] for i in range(N)]
        for a,r in enumerate(s.reps):
            La=[[s.so(a,(u-s.th[a])*(i==j)*I+h*r(i,j)) for j in range(N)] for i in range(N)]
            T=[[sum(T[i][k]@La[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
        s.cache[u]=T; return T
    def t(s,a,u):
        D=N**s.L
        if a==0: return np.eye(D)
        out=np.zeros((D,D))
        for I_ in itertools.combinations(range(N),a):
            for J in itertools.combinations(range(N),a):
                dg=np.linalg.det(G[np.ix_(J,I_)])
                if abs(dg)<1e-14: continue
                q=np.zeros((D,D))
                for p in itertools.permutations(range(a)):
                    sg=np.linalg.det(np.eye(a)[list(p)])
                    q+=sg*functools.reduce(np.matmul,[s.T(u-k*h)[I_[p[k]]][J[k]] for k in range(a)])
                out+=q*dg
        return out
th=[0.2857,-0.6,0.731,-1.37]
C4=Chain([rA,rF,rA,rF],th); C2=Chain([rA,rF],th[:2])
rng=np.random.default_rng(5); X=0.3*rng.standard_normal((N,N)); w,V=np.linalg.eig(X); g=(V@np.diag(np.exp(w))@np.linalg.inv(V)).real
gA=g; gF=np.linalg.inv(g).T
e1=I[0]; e4=I[3]
Om2=np.kron(e1@gA,e4@gF); Om4=np.kron(Om2,Om2)
tA4={a:C4.t(a,th[0]-h) for a in range(4)}; tF4=C4.t(1,th[1])
tA2={a:C2.t(a,th[0]-h) for a in range(4)}; tF2=C2.t(1,th[1])
for a in range(4):
    for k in range(4):
        v=Om4@tA4[a]@np.linalg.matrix_power(tF4,k)
        M=v.reshape(16,16); s=np.linalg.svd(M,compute_uv=False)
        target=np.kron(Om2@tA2[a]@np.linalg.matrix_power(tF2,k),Om2)
        c=np.vdot(target,v)/np.vdot(target,target); r=np.linalg.norm(v-c*target)/np.linalg.norm(v)
        print("a=%d k=%d  rank1-ratio %.1e  vs L2(x)Om: resid %.1e  c=%.6g"%(a,k,s[1]/s[0],r,c))
