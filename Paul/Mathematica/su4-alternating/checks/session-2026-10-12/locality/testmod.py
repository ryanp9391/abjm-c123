from modlib import *
import time
p = PRIMES[0]; print(PRIMES)
L = 4
mc = ModChain(L, TH6, H, Z4, p)
ch = build_chain(L, dt=Fraction)  # exact
# a = 1/3 W state exact
a = Fraction(1,3)
def Wex(L,a,a1=None):
    out = np.array([Fraction(1)],dtype=object)
    for pp in range(L//2):
        aa = a1 if (pp==0 and a1 is not None) else a
        v = np.array([Fraction(1),0,0,aa],dtype=object); u=np.array([-aa,0,0,Fraction(1)],dtype=object)
        out = np.kron(np.kron(out,v),u)
    return out
psi = Wex(L,a)
u = TH6[0] + Fraction(1,2)*(-2)
for which in (1,2,3):
    uu = TH6[1] + Fraction(1,2)*3
    ex = ch.apply(which, uu, psi)
    # mod
    Psim = np.array([[fm(x,p)] for x in psi], dtype=float)
    mm = mc.apply(which, uu, Psim)
    exm = np.array([fm(x,p) for x in ex], dtype=float)
    print(which, np.array_equal(exm, mm[:,0]), ex[:3])
# float comparison with exact
chf = build_chain(L)
psif = np.array([float(x) for x in psi])
for which in (1,2,3):
    uu = TH6[1]+Fraction(1,2)*3
    print(which, np.abs(np.array([float(x) for x in ch.apply(which,uu,psi)]) - chf.apply(which,uu,psif)).max())
