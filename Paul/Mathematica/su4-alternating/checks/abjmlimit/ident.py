from words import *
F = ModP(); h = H; p = F.p
tF = Fr(-3,5); tA = tF + 2*h
a, b = Fr(-31,40), Fr(2,9)
for name, th in [('only 1-2 resonant', [tA, tF, tA+a, tF+b]), ('pair-by-pair', [tA, tF, tA+a, tF+a]), ('pair-periodic', [tA,tF,tA,tF])]:
    ch = Chain(F, TYPES, th, h, G)
    u = th[0]-h
    t1 = ch.t(1,u); t2 = ch.t(2,u); t3 = ch.t(3,u); Y = ch.t(1, th[1])
    print(name)
    print('  t3(thA-h)==0 :', np.count_nonzero(t3)==0)
    print('  t2(thA-h) == t1(thA-h) t1(thF) :', np.count_nonzero((t2 - F.mm(t1,Y))%p)==0)
    print('  t2(thA-h) == t1(thF) t1(thA-h) :', np.count_nonzero((t2 - F.mm(Y,t1))%p)==0)
    vac = vacuum(TYPES, F)
    print('  on <0| only:', np.count_nonzero((vac@t2 - vac@F.mm(t1,Y))%p)==0)
