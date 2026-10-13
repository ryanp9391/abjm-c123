import numpy as np
d = np.load('states_std.npz')
WA1, WA3, WF2 = d['WA1'], d['WA3'], d['WF2']
for name, WA in (('site1', WA1), ('site3', WA3)):
    t1, t2, t3 = WA[1], WA[2], WA[3]
    r2 = t2 / t1 ** 2; r3 = t3 / t1 ** 3; r32 = t3 * t1 / t2 ** 2
    print(name, 'tau2/tau1^2 spread', np.std(r2) / np.abs(np.mean(r2)), ' tau3/tau1^3', np.std(r3) / np.abs(np.mean(r3)), 'tau3 tau1/tau2^2', np.std(r32)/np.abs(np.mean(r32)))
    # linear relation test: rank of [1, t1, t2, t3, t1^2, t1^3, t1 t2] across states
    cols = np.array([np.ones_like(t1), t1, t2, t3, t1**2, t1**3, t1*t2, t2**2]).T
    s = np.linalg.svd(cols / np.abs(cols).max(axis=0), compute_uv=False)
    print('  singular values', s / s[0])
    print('  # distinct tau1 values', len(np.unique(np.round(t1, 6))), ' tau2', len(np.unique(np.round(t2, 6))), ' tau3', len(np.unique(np.round(t3, 6))))
