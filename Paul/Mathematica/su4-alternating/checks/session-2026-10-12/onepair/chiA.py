from common import *
def Qt(s, j, x): return z[j]**(x / h) * qj(s, j, x)
xs = [th1 - b * h for b in range(4)]
for s in states:
    Qrows = [[Qt(s, i, x) for x in xs] for i in range(3)]
    D = [det([[row[b] for b in range(4) if b != a] for row in Qrows]) for a in range(4)]
    tt = [D[a] / D[0] for a in range(4)]          # tau-tilde_a
    chi = [(-1)**b * tt[b] for b in range(4)]      # coefficients of y^3, y^2, y, 1
    rts = mp.polyroots(chi, maxsteps=200, extraprec=200)
    R = tau(s, 1, th2) / (h * th)
    zs = [mp.nstr(x, 8) for x in z]
    print(s['label'], 'R =', mp.nstr(R, 8), ' roots of chi_A:', [mp.nstr(x, 8) for x in rts], ' tau1(th1)/(h th) =', mp.nstr(tau(s, 1, th1) / (h * th), 8))
print('z =', [mp.nstr(x, 8) for x in z])
