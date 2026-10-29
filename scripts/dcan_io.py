"""Read DCAN motion masks (MATLAB v5 .mat) with the standard library only."""
import struct
import zlib

TYPES = {1: ('b', 1), 2: ('B', 1), 3: ('h', 2), 4: ('H', 2), 5: ('i', 4),
         6: ('I', 4), 7: ('f', 4), 9: ('d', 8), 12: ('q', 8), 13: ('Q', 8)}


def elems(buf):
    i = 0
    while i < len(buf):
        t, nb = struct.unpack('<II', buf[i:i + 8])
        if t >> 16:
            nb = t >> 16
            t = t & 0xffff
            body = buf[i + 4:i + 4 + nb]
            i += 8
        else:
            body = buf[i + 8:i + 8 + nb]
            i += 8 + ((nb + 7) // 8) * 8
        if t == 15:
            for x in elems(zlib.decompress(body)):
                yield x
        else:
            yield (t, body)


def parse_matrix(body):
    subs = list(elems(body))
    cls = subs[0][1][0]
    dims = struct.unpack('<%di' % (len(subs[1][1]) // 4), subs[1][1])
    if cls == 1:  # cell array
        return ('cell', dims, [parse_matrix(b) for t, b in subs[3:] if t == 14])
    if cls == 2:  # struct
        nlen = struct.unpack('<i', subs[3][1][:4])[0]
        names = subs[4][1]
        fn = [names[i:i + nlen].rstrip(b'\x00').decode()
              for i in range(0, len(names), nlen)]
        return ('struct', dims, {nm: parse_matrix(subs[5 + i][1])
                                 for i, nm in enumerate(fn)})
    dt, raw = subs[3]
    if dt in TYPES:
        fmt, sz = TYPES[dt]
        n = len(raw) // sz
        vals = list(struct.unpack('<%d%s' % (n, fmt), raw[:n * sz]))
    elif dt == 2:
        vals = raw.decode(errors='ignore')
    else:
        vals = raw
    return ('num', dims, vals)


def load_mask(path):
    data = open(path, 'rb').read()
    for t, body in elems(data[128:]):
        if t == 14:
            subs = list(elems(body))
            if subs[2][1].decode(errors='ignore') == 'motion_data':
                return parse_matrix(body)
    raise ValueError('motion_data not found in ' + path)
