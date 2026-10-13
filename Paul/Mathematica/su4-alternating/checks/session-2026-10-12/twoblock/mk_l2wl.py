import numpy as np, zlib
import wxf
raw = open(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf', 'rb').read()
print('compressed', len(raw), 'decompressed', len(zlib.decompress(raw[3:])))
D = wxf.load(r'C:/Users/pryan/git-projects/abjm-c123/Paul/Mathematica/su4-alternating/data/abjm-hw-states.wxf')
def num(x):
    x = complex(x)
    re, im = repr(x.real), repr(x.imag)
    f = lambda s: s.replace('e', '*^')
    return f(re) if x.imag == 0 else f'Complex[{f(re)}, {f(im)}]'
out = []
for s in D['L2']:
    sa = s['psi']; cols = np.array(sa.args[3][1][1]).ravel(); vals = sa.args[3][2]
    rules = ', '.join(f'{int(c)} -> {num(v)}' for c, v in zip(cols, vals))
    C = s['C']
    out.append('<|"w" -> {%s}, "dyn" -> {%s}, "E" -> %s, "C" -> {%s}, "psi" -> SparseArray[{%s}, 256]|>' % (
        ', '.join(str(int(x)) for x in s['w']), ', '.join(str(int(x)) for x in s['dyn']), num(s['E']), ', '.join(num(c) for c in C), rules))
txt = 'tbL2 = {' + ',\n  '.join(out) + '};'
open('l2_states.wl', 'w').write(txt)
print(len(txt)); print(txt[:600])
