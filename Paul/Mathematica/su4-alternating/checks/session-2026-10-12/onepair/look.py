from common import *
for s in states[:3] + states[12:13]:
    print(s['label'])
    for j in range(4):
        cs = [M(c) for c in s['q'][j]]
        deg = max(m for m in range(4) if abs(cs[m]) > 1e-30)
        rts = mp.polyroots(cs[:deg + 1][::-1]) if deg > 0 else []
        print('  q%d deg %d roots' % (j + 1, deg), [mp.nstr(x, 12) for x in rts], ' (roots - th1)/h:', [mp.nstr((x - th1) / h, 10) for x in rts], ' (roots-th2)/h:', [mp.nstr((x - th2) / h, 10) for x in rts])
