import tempfile
import sys; sys.path.insert(0,'.')
import numpy as np, wxf, itertools, pickle, os
PATH=r'C:\Users\pryan\git-projects\abjm-c123\Paul\Mathematica\su4-alternating\data\abjm-hw-states.wxf'
def dense(psi, N):
    # psi: ('SparseArray', ['Automatic', dims, bg, [1, [rowptr, cols], vals]]) or ndarray
    if isinstance(psi,np.ndarray): return np.asarray(psi,dtype=complex).reshape(-1)
    h,args=psi
    assert h=='SparseArray'
    dims=int(args[1][0]); bg=args[2]; rp,cols=args[3][1]; vals=args[3][2]
    v=np.full(dims,bg,dtype=complex)
    cols=np.asarray(cols).reshape(-1)-1
    v[cols]=np.asarray(vals,dtype=complex)
    return v
def load():
    cache=os.path.join(tempfile.gettempdir(),'twoblock_states.pkl')
    if os.path.exists(cache): return pickle.load(open(cache,'rb'))
    d=wxf.parse(open(PATH,'rb').read())
    out={}
    for L in (2,3,4):
        lst=[]
        for s in d['L%d'%L]:
            t=dict(s); t['psi']=dense(s['psi'],4**(2*L)); t['w']=np.array(s['w']); t['dyn']=np.array(s['dyn']); lst.append(t)
        out[L]=lst
    pickle.dump(out,open(cache,'wb'))
    return out
if __name__=='__main__':
    S=load()
    for L in (2,3,4):
        print(L,len(S[L]))
        for k,s in enumerate(S[L][:8]):
            print(k,s['w'],s['dyn'],round(s['E'],6),np.count_nonzero(np.abs(s['psi'])>1e-14),np.linalg.norm(s['psi']),[complex(np.round(c,4)) for c in s['C']])
