# Fallback python exploration for twoblock-polarisation-forms.wb (written while the Wolfram kernel was down).
# Run order (from this folder, plain `python <script>`): states.py (WXF -> pickle in the temp dir), run_all.py (fits),
# polyverify.py 2 3 4, invdim.py 2 3 4, invdim_all.py, formula.py, symm.py, gen_json.py (writes data/abjm-twoblock-forms.json), sympycheck.py, nbdata.py.
import os, tempfile
import pickle, time, sys
from inv import *
res={}
t0=time.time()
for L in (2,3,4):
    for k,s in enumerate(S[L]):
        for M in range(L+1):
            for bps in (True,False):
                if L==4 and not bps and M not in (0,1,2,3,4): continue
                r=fit(L,s,M,bps)
                # store compact
                res[(L,k,M,bps)]=dict(nm=r['nm'],rank=r['rank'],Fmax=r['Fmax'],res=r['res'] if r['nm'] else None,
                    c=(r['c'] if r['nm'] else None), mons=(r['mons'] if r['nm'] else None))
    print(L,'done',time.time()-t0,flush=True)
    pickle.dump(res,open(os.path.join(tempfile.gettempdir(),'twoblock_fits.pkl'),'wb'))
