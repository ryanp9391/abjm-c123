import struct, zlib, numpy as np
class R:
    def __init__(s,b): s.b=b; s.p=0
    def byte(s): v=s.b[s.p]; s.p+=1; return v
    def varint(s):
        r=0; sh=0
        while True:
            c=s.byte(); r|=(c&127)<<sh; sh+=7
            if c<128: return r
    def take(s,n): v=s.b[s.p:s.p+n]; s.p+=n; return v
class Sym(str): pass
def parse(data):
    assert data[:2]==b'8:' or data[:3]==b'8C:', data[:5]
    if data[:3]==b'8C:': body=zlib.decompress(data[3:])
    else: body=data[2:]
    r=R(body)
    return rd(r)
NT={0:('b',1),1:('<h',2),2:('<i',4),3:('<q',8),16:('B',1),17:('<H',2),18:('<I',4),19:('<Q',8),34:('<f',4),35:('<d',8)}
def rd(r):
    t=r.byte()
    c=chr(t)
    if c=='f':
        n=r.varint(); h=rd(r); args=[rd(r) for _ in range(n)]
        hn=str(h)
        if hn=='List': return args
        if hn=='Rational': 
            from fractions import Fraction
            return Fraction(int(args[0]),int(args[1]))
        if hn=='Complex': return complex(float(args[0]),float(args[1]))
        return (hn,args)
    if c=='s': n=r.varint(); return Sym(r.take(n).decode())
    if c=='S': n=r.varint(); return r.take(n).decode('utf8')
    if c in 'Ci': return struct.unpack('<b',r.take(1))[0]
    if c=='j': return struct.unpack('<h',r.take(2))[0]
    if c=='k': return struct.unpack('<i',r.take(4))[0]
    if c=='L': return struct.unpack('<q',r.take(8))[0]
    if c=='r': return struct.unpack('<d',r.take(8))[0]
    if c=='I': n=r.varint(); return int(r.take(n).decode())
    if c=='R': n=r.varint(); return float(r.take(n).decode().split('`')[0])
    if c=='A':
        n=r.varint(); d={}
        for _ in range(n):
            k=r.byte(); assert chr(k) in '-:'
            key=rd(r); val=rd(r); d[key]=val
        return d
    if t==193:
        ty=r.byte(); rank=r.varint(); dims=[r.varint() for _ in range(rank)]
        N=int(np.prod(dims))
        if ty in NT:
            f,sz=NT[ty]; a=np.frombuffer(r.take(N*sz),dtype=np.dtype(f)).reshape(dims); return a
        if ty==51: a=np.frombuffer(r.take(N*8),dtype='<f4').reshape(dims+[2]); return a[...,0]+1j*a[...,1]
        if ty==52: a=np.frombuffer(r.take(N*16),dtype='<f8').reshape(dims+[2]); return a[...,0]+1j*a[...,1]
        raise ValueError(ty)
    if t==194:
        ty=r.byte(); rank=r.varint(); dims=[r.varint() for _ in range(rank)]
        N=int(np.prod(dims))
        if ty in NT:
            f,sz=NT[ty]; return np.frombuffer(r.take(N*sz),dtype=np.dtype(f)).reshape(dims)
        if ty==51: a=np.frombuffer(r.take(N*8),dtype='<f4').reshape(dims+[2]); return a[...,0]+1j*a[...,1]
        if ty==52: a=np.frombuffer(r.take(N*16),dtype='<f8').reshape(dims+[2]); return a[...,0]+1j*a[...,1]
        raise ValueError(ty)
    raise ValueError('token %d %s at %d'%(t,c,r.p))
