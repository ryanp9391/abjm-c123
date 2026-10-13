import struct, zlib, numpy as np
from fractions import Fraction
def parse(data):
    assert data[:2]==b'8:' or data[:3]==b'8C:', data[:5]
    if data[:3]==b'8C:': data=zlib.decompress(data[3:])
    else: data=data[2:]
    pos=[0]
    def varint():
        r=0; sh=0
        while True:
            b=data[pos[0]]; pos[0]+=1
            r|=(b&0x7f)<<sh; sh+=7
            if not b&0x80: return r
    def rd(n):
        s=data[pos[0]:pos[0]+n]; pos[0]+=n; return s
    types={0:('b',1),1:('h',2),2:('i',4),3:('q',8),4:('B',1),5:('H',2),6:('I',4),7:('Q',8),16:('f',4),17:('d',8),51:('f',4),52:('d',8)}
    def sparse(args):
        # SparseArray[Automatic, dims, default, {1, {ptrs, cols}, vals}]
        dims=[int(x) for x in np.atleast_1d(args[1])]; dflt=args[2]; body=args[3]
        assert body[0]==1
        ptrs=np.array(body[1][0]).astype(int); cols=np.array(body[1][1]).astype(int)-1; vals=np.array(body[2])
        if len(dims)==1:
            out=np.full(dims[0],dflt,dtype=vals.dtype if vals.size else float); 
            # rank-1: ptrs=[0,n], cols shape (n,1)
            out[cols.reshape(-1)]=vals; return out
        out=np.full(dims,dflt,dtype=vals.dtype if vals.size else float)
        for r in range(dims[0]):
            for k in range(ptrs[r],ptrs[r+1]): out[(r,)+tuple(np.atleast_1d(cols[k]))]=vals[k]
        return out
    def tok():
        t=rd(1)[0]
        if t==102:
            n=varint(); head=tok(); args=[tok() for _ in range(n)]
            if head=='Rational': return Fraction(args[0],args[1])
            if head=='Complex': return complex(args[0],args[1])
            if head=='List': return args
            if head=='SparseArray': return sparse(args)
            return (head,args)
        if t==115: n=varint(); return rd(n).decode()
        if t==83: n=varint(); return rd(n).decode('utf8')
        if t in (67,105): return struct.unpack('<b',rd(1))[0]
        if t==106: return struct.unpack('<h',rd(2))[0]
        if t==107: return struct.unpack('<i',rd(4))[0]
        if t==76: return struct.unpack('<q',rd(8))[0]
        if t==114: return struct.unpack('<d',rd(8))[0]
        if t==73: n=varint(); return int(rd(n).decode())
        if t==82: n=varint(); return float(rd(n).decode().split('`')[0])
        if t==66: n=varint(); return rd(n)
        if t==65:
            n=varint(); d={}
            for _ in range(n):
                k=rd(1)[0]; key=tok(); val=tok(); d[key]=val
            return d
        if t in (193,194):
            ty=rd(1)[0]; rank=varint(); dims=[varint() for _ in range(rank)]
            f,sz=types[ty]; cnt=int(np.prod(dims)) if dims else 1
            mult=2 if ty>=51 else 1
            arr=np.array(struct.unpack('<%d%s'%(cnt*mult,f),rd(cnt*mult*sz)))
            if mult==2: arr=arr[0::2]+1j*arr[1::2]
            return arr.reshape(dims)
        raise ValueError('token %d at %d'%(t,pos[0]))
    return tok()
def load(path): return parse(open(path,'rb').read())
if __name__=='__main__':
    import sys
    d=load(sys.argv[1])
    def show(x,ind=0,maxd=3):
        if isinstance(x,dict):
            for k,v in x.items():
                print(' '*ind+repr(k)+':',type(v).__name__, (v.shape if isinstance(v,np.ndarray) else ''))
                if ind<6: show(v,ind+2)
        elif isinstance(x,list):
            print(' '*ind+'list',len(x))
            if x and ind<6: show(x[0],ind+2)
        elif isinstance(x,str): print(' '*ind+x[:2500])
    show(d)
