# sanity: data psi are eigenvectors of the untwisted homogeneous chain (ket action), and C values reproduce
import numpy as np
import untw, eng
ch, z = untw.make_chain(0.0, [0, 0, 0, 0], [0, 0, 0])
for i in (0, 16, 18, 20, 22):
    sec = untw.SECT[tuple(untw.W0[i])]
    t, r = untw.eigval(ch, untw.PSI0[i], sec, 1, 0.37)
    t3, r3 = untw.eigval(ch, untw.PSI0[i], sec, 3, -0.29, True)
    # C_2 = <vac W-1^2|psi>, C_4 = <vac^2 W-1|psi>
    vac = eng.VAC
    B2 = eng.prod_state(vac + eng.W(-1) * 2, 'c').reshape(-1); B4 = eng.prod_state(vac * 2 + eng.W(-1), 'c').reshape(-1)
    print(i, 'res t1 %.1e t3 %.1e' % (r, r3), ' C2: %s vs data %s' % (np.round(B2 @ untw.PSI0[i], 8), np.round(untw.C0[i][2], 8)),
          ' C4: %s vs %s' % (np.round(B4 @ untw.PSI0[i], 8), np.round(untw.C0[i][4], 8)))
