import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from symR6 import *
u3 = sp.Symbol('u3')
solve([1, 1, 0, 0], [0, 0, u3, u4], 'R2 (Levi rep): v = (1,1,0,0), u = (0,0,u3,u4)')
