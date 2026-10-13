from words import *
F = ModP(); h = H
def rk(th):
    ch = Chain(F, TYPES, th, h, G)
    return rank_mod(words(ch))
tF = Fr(-3,5); tA = tF + 2*h
a, b = Fr(-31,40), Fr(2,9)
cases = [
 ('generic d3,d4 (only th1-th2=2h)', [tA, tF, tA+a, tF+b]),
 ('d3=d4 (pair-by-pair: th1-th2=th3-th4=2h)', [tA, tF, tA+a, tF+a]),
 ('d4=0 (th1-th2=th1-th4=2h; th3 generic)', [tA, tF, tA+a, tF]),
 ('d3=0 (th1=th3; th1-th2=th3-th2=2h)', [tA, tF, tA, tF+b]),
 ('d3=d4=0 (pair-periodic, ABJM)', [tA, tF, tA, tF]),
 ('only th1-th4=2h', [tA, tF+b, tA+a, tF]),
 ('only th3-th2=2h', [tA+b, tF, tA+a+b, tF+a+2*b] if False else [tF+b, tF, tF+2*h, tF+a]),
 ('only th3-th4=2h', [tA+b, tF, tA+a, tF+a]),
 ('th1-th2 = 2h and th3-th4 = 2h and th1-th4=2h? (needs d3=d4=0)', [tA, tF, tA, tF]),
 ('pair-periodic generic d', [tF+Fr(1,3), tF, tF+Fr(1,3), tF]),
 ('th1=th3 only (generic d)', [tF+Fr(1,3), tF, tF+Fr(1,3), tF+b]),
 ('th2=th4 only (generic d)', [tF+Fr(1,3), tF, tF+Fr(1,3)+a, tF]),
]
for name, th in cases:
    print(name, [str(x) for x in th], 'rank', rk(th), flush=True)
