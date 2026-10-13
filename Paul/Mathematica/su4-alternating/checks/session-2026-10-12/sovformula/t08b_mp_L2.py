"""High-precision L = 2 check of the A-site (odd-bridge) word on <0| and of the norm identity |<0|Psi>|^2 = |lambda^w N|."""
import numpy as np, sys, json
from mpu import *
from untw import load_data_states
DATA = load_data_states(2)
phis = {'phiB': ['1/5', '-7/10', '1/2', '13/10']}
thdir = ['1/3', '-4/5']
y0 = Fr(-13, 10)
for pname, phi in phis.items():
    for eps in (Fr(1, 10 ** 6),):
        C = Ctx(2, eps, [Fr(p) for p in phi], [Fr(t) for t in thdir], dps=80)
        lam = C.lam; z = C.z; c4 = z[0] * z[1] * z[2] * z[3]
        th = C.th; uF, u0 = th[1], th[1] + H
        l3 = Fr(1)
        for s in (2, 3):
            t = th[s]
            fac = [u0 - t - j * H for j in range(3)] if s % 2 == 0 else [u0 - t - 3 * H] + [u0 - t - j * H for j in range(2)]
            for f in fac:
                if f != 0: l3 *= f
        alpha = c4 / z[0] + 1 / (lam[0] * y0); beta = -1 / (12 * H ** 5 * l3)
        out = []
        for st in DATA[:3]:
            w = st['w']
            T3 = C.secmat(w, 3, u0, reg=True)
            Mg = [[a + Fr(61, 100) * b for a, b in zip(r1, r2)] for r1, r2 in zip(C.secmat(w, 1, uF), C.secmat(w, 3, Fr(-29, 100)))]
            M2 = C.secmat(w, 2, Fr(83, 100))
            Mg = [[a + Fr(23, 100) * b for a, b in zip(r1, r2)] for r1, r2 in zip(Mg, M2)]
            psi, psiL, ov, gap = C.eig(w, Mg, st['psi'])
            tv = matvec(T3, psi); tau3 = dot([mp.conj(x) for x in psi], tv)
            bra0 = C.prodbra([lam, [1] * 4] * 2, w); ket0b = C.prodbra([[1] * 4, lam] * 2, w)
            o0 = dot(bra0, psi); o0b = dot(psiL, ket0b); nn = dot(psiL, psi)
            Nsov = o0 * o0b / nn
            lw = mp.fprod([tomp(lam[j]) ** w[j] for j in range(4)])
            P = tomp(alpha) + tomp(beta) * tau3
            out.append(dict(st=st, o0=o0, Nsov=Nsov, lw=lw, P=P, herm=abs(o0) ** 2 / abs(lw * Nsov), ov=ov))
        vac = out[0]
        print(f'{pname} eps={float(eps):.0e}')
        for r in out:
            Dw = (tomp(lam[1] / lam[0])) ** r['st']['w'][1] * (tomp(lam[2] / lam[0])) ** r['st']['w'][2]
            Cdir = r['P'] * r['o0'] / (vac['P'] * vac['o0']) / Dw
            Cform = (r['P'] / vac['P']) * mp.sqrt(r['lw'] / vac['lw'] * r['Nsov'] / vac['Nsov']) / Dw
            print(f"   {list(r['st']['dyn'])} herm={mp.nstr(r['herm'], 12)}  C_direct={mp.nstr(Cdir, 10)} |C_direct|={mp.nstr(abs(Cdir), 10)}  |C_SoVformula|={mp.nstr(abs(Cform), 10)}  <0|Psi>/lam1^2={mp.nstr(r['o0'] / tomp(lam[0]) ** 2, 6)}")
        sys.stdout.flush()
