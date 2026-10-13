import numpy as np
exec(open('s04_Cop.py').read().split('results = []')[0])
C1, C2 = C_cand(u1, False, 1, 0, True), C_cand(u2, False, 1, 0, True)
print('[C,B] rel:', np.abs(C1 @ Bu - Bu @ C1).max() / (np.abs(C1).max() * np.abs(Bu).max()))
# eigen-decomposition of C at u1 (right eigenvectors = candidate right SoV basis)
ev, Yr = np.linalg.eig(C1)
# check that C(u) eigenvalues factorise: compute eigenvalue of C(u) at several u for each eigvec, fit polynomial roots
us = [Fr(k, 7) for k in range(-12, 13)]
vals = []
for u in us:
    Cu = C_cand(u, False, 1, 0, True)
    vals.append(np.einsum('ij,ij->j', np.linalg.pinv(Yr).T, Cu @ Yr))
vals = np.array(vals)
deg = None
uf = np.array([float(u) for u in us])
for d in range(8, 26):
    coef = np.polyfit(uf, vals[:, 0].real, d)
    if np.abs(np.polyval(coef, uf) - vals[:, 0].real).max() < 1e-6 * np.abs(vals[:, 0]).max(): deg = d; break
print('degree of C eigenvalue polynomial (state 0):', deg)
th = [float(x) for x in th0]
for n in range(3):
    coef = np.polyfit(uf, vals[:, n].real, deg)
    print('roots (minus thetas):', np.round(np.sort_complex(np.roots(coef)), 4))
print('thetas', th, 'h', h)
# left B-eigen words vs right C-eigenvectors: overlap matrix <x|y>
M, vac = sov_words(Chain(R))
N = M @ Yr
print('nonzeros of <x|y> (|.|>1e-8 max):', int((np.abs(N) > 1e-8 * np.abs(N).max()).sum()))
np.save('Cright.npy', Yr)
