"""load the untwisted HW cyclic states at L = 3 pairs from data/abjm-hw-states.wxf (key 'L3'), cached as npz"""
import numpy as np, os
DATA = "C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf"
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hw_L3.npz')
def load():
    if os.path.exists(CACHE):
        d = np.load(CACHE, allow_pickle=True); return d['psi'], d['w'], d['dyn'], d['E'], d['C']
    import wxf
    d = wxf.load(DATA)['L3']
    psi = np.array([np.array(s['psi'], dtype=complex) for s in d])
    w = np.array([np.array(s['w'], dtype=int) for s in d]); dyn = np.array([np.array(s['dyn'], dtype=int) for s in d])
    E = np.array([float(s['E']) for s in d]); C = np.array([np.array(s['C'], dtype=complex) for s in d])
    np.savez(CACHE, psi=psi, w=w, dyn=dyn, E=E, C=C)
    return psi, w, dyn, E, C
