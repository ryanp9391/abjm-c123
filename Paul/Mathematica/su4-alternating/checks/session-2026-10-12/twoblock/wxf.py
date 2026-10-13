import zlib, struct
import numpy as np
from fractions import Fraction

class Sym(str): pass
class Fn:
    def __init__(s, head, args): s.head, s.args = head, args
    def __repr__(s): return f"{s.head}[{', '.join(map(repr, s.args[:5]))}{'...' if len(s.args)>5 else ''}]"

def varint(b, i):
    r = 0; sh = 0
    while True:
        c = b[i]; i += 1
        r |= (c & 0x7f) << sh; sh += 7
        if c < 128: return r, i

NA = {0x00:'i1',0x01:'<i2',0x02:'<i4',0x03:'<i8',0x10:'u1',0x11:'<u2',0x12:'<u4',0x13:'<u8',0x22:'<f4',0x23:'<f8',0x33:'<c8',0x34:'<c16'}

def parse(b, i):
    t = b[i]; i += 1
    c = chr(t)
    if c == 'f':
        n, i = varint(b, i); head, i = parse(b, i)
        args = []
        for _ in range(n):
            a, i = parse(b, i); args.append(a)
        return post(Fn(head, args)), i
    if c in 'sSIR':
        n, i = varint(b, i); s = b[i:i+n].decode('utf8'); i += n
        if c == 's': return Sym(s), i
        if c == 'S': return s, i
        if c == 'I': return int(s), i
        if c == 'R':
            return float(s.split('`')[0].replace('*^','e')), i
    if c == 'C': return struct.unpack('<b', b[i:i+1])[0], i+1
    if c == 'j': return struct.unpack('<h', b[i:i+2])[0], i+2
    if c == 'i': return struct.unpack('<i', b[i:i+4])[0], i+4
    if c == 'L': return struct.unpack('<q', b[i:i+8])[0], i+8
    if c == 'r': return struct.unpack('<d', b[i:i+8])[0], i+8
    if c == 'A':
        n, i = varint(b, i); d = {}
        for _ in range(n):
            i += 1  # '-' or ':'
            k, i = parse(b, i); v, i = parse(b, i); d[k] = v
        return d, i
    if t in (0xC1, 0xC2):
        typ = b[i]; i += 1
        r, i = varint(b, i); dims = []
        for _ in range(r):
            d, i = varint(b, i); dims.append(d)
        dt = np.dtype(NA[typ]); cnt = int(np.prod(dims))
        arr = np.frombuffer(b[i:i+cnt*dt.itemsize], dtype=dt).reshape(dims); i += cnt*dt.itemsize
        return arr, i
    raise ValueError(f"unknown token {t} at {i}")

def post(f):
    h = f.head
    if h == 'Complex' and len(f.args) == 2:
        return complex(f.args[0], f.args[1]) if not all(isinstance(a,int) for a in f.args) else complex(*f.args)
    if h == 'Rational': return Fraction(f.args[0], f.args[1])
    if h == 'List': return f.args
    return f

def load(path):
    b = open(path,'rb').read()
    if b[:3] == b'8C:': b = b'8:' + zlib.decompress(b[3:])
    assert b[:2] == b'8:'
    v, i = parse(b, 2)
    return v
