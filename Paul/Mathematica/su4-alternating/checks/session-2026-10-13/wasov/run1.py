from wa_sov_core import *
ch, W, Y, Yk, G, M = build()
print('Gram nnz', sum(1 for r in G for x in r if x != 0), 'M nnz', len(M))
