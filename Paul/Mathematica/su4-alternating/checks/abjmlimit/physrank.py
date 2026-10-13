from words import *
F = ModP(); h = H; p = F.p
tF = Fr(-3,5); tA = tF + 2*h
th = [tA, tF, tA, tF]
ch = Chain(F, TYPES, th, h, G)
def wshift(e1, z2, e3, z4):
    vac = vacuum(TYPES, F)
    A1 = [ch.t(a, tA-h+e1) for a in range(4)]; A3 = [ch.t(a, tA-h+e3) for a in range(4)]
    T2 = ch.t(1, tF+z2); T4 = ch.t(1, tF+z4)
    rows = []
    for (a1,k2,a3,k4) in WIDX:
        v = vac @ A1[a1] % p
        for _ in range(k2): v = v @ T2 % p
        v = v @ A3[a3] % p
        for _ in range(k4): v = v @ T4 % p
        rows.append(v)
    return np.array(rows)
for sh in [(Fr(1,7),Fr(2,9),Fr(-1,5),Fr(3,11)), (Fr(0),Fr(0),Fr(1,7),Fr(2,9)), (Fr(1,7),Fr(0),Fr(-1,5),Fr(0)), (Fr(0),Fr(2,9),Fr(0),Fr(3,11))]:
    print('physical point, shifts (eta1, zeta2, eta3, zeta4) =', [str(x) for x in sh], ' rank', rank_mod(wshift(*sh)), flush=True)
    ch._tcache.clear(); ch._Tcache.clear()
