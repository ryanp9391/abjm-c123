import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'g2test.py')).read().split('for a in range(4):\n    for k')[0])
inv=np.linalg.inv
ops4={}; ops2={}
for name,(a,u) in {'t1(th1)':(1,th[0]),'t3(th1-h)':(3,th[0]-h),'t1(th2)':(1,th[1]),'t3(th2+3h)':(3,th[1]+3*h),
                   't2(th1)':(2,th[0]),'t2(th1-h)':(2,th[0]-h),'t2(th2+h)':(2,th[1]+h),'t2(th2+2h)':(2,th[1]+2*h),'t2(th2)':(2,th[1]),'t2(th1+h)':(2,th[0]+h)}.items():
    ops4[name]=C4.t(a,u); ops2[name]=C2.t(a,u)
for k in list(ops4):
    if k in ('t1(th1)','t1(th2)'):
        ops4[k+'^-1']=inv(ops4[k]); ops2[k+'^-1']=inv(ops2[k])
def loc(v4,v2):
    target=np.kron(v2,Om2); c=np.vdot(target,v4)/np.vdot(target,target)
    return np.linalg.norm(v4-c*target)/np.linalg.norm(v4)
names=list(ops4)
print("singles:")
for n in names: print("  %-14s %.1e"%(n,loc(Om4@ops4[n],Om2@ops2[n])))
print("pairs (localising only):")
for n1 in names:
    for n2 in names:
        if n1>=n2: continue
        r=loc(Om4@ops4[n1]@ops4[n2],Om2@ops2[n1]@ops2[n2])
        if r<1e-9: print("  %s * %s"%(n1,n2))
