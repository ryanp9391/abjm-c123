import os; D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,'scan_classify.py')).read().split('want=[')[0])
k1,k2=0.37,-0.81; v=e1@(I-k1*G.T)@(I-k2*G.T); refs=[v,e4,v,e4]; Om=prod(refs)
a=Om@ops['t3(th1-1 h)']; b=Om@ops['t1(th1+0 h)^-1']
print("t3(th1-h) vs t1(th1)^-1 on Omega_2 proportional:", abs(abs(np.vdot(a,b))-np.linalg.norm(a)*np.linalg.norm(b))<1e-9*np.linalg.norm(a)*np.linalg.norm(b))
x=prod([v@np.linalg.inv(G).T,e4,v,e4]); print("  = rho1(G^-T) on site 1:", np.linalg.matrix_rank(np.array([b,x]),1e-9)==1)
sub=[n for n in names if any(n.startswith(p) for p in ('t1(th','t3(th'))]
sub=[n for n in sub if ('+0 h' in n and n.startswith('t1')) or ('-1 h' in n and n.startswith('t3'))]
print("subset:",sub)
import itertools as it
found=[]
for w in it.combinations_with_replacement(sub,3):
    psi=Om
    for n in w: psi=psi@ops[n]
    cl=classify(refs,psi,Om)
    if cl in('zero','prop'): continue
    for B in loc_blocks(refs,psi):
        if len(B)>=2: found.append((B,w,cl))
for B in [(0,1,2),(1,2,3),(0,1,3),(0,2,3),(0,1),(2,3),(0,3),(1,2),(0,2),(1,3)]:
    L=[(w,c) for b,w,c in found if b==B]; P=[x for x in L if x[1]=='product']
    print("sites",[s+1 for s in B],": %d triple words, %d product"%(len(L),len(P)), P[:4])
