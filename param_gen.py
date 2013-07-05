import math
import constants
import operator
from bitarray import bitarray

def tuples(encoder, X):
    '''
    Generate encoding tuples as described here http://tools.ietf.org/html/rfc6330#section-5.3.5.4
    '''
    A = 53591 + encoder.J*997
    if A%2 == 0:
        A += 1
    B = 10267*(encoder.J+1)
    y = int((B + X*A) % math.pow(2,32))
    v = random(y, 0, math.pow(2, 20))
    d = degree(v, encoder.W)
    if d is None:
        raise Exception("Error producing value d from y=%d and W=%d" % (y, encoder.W))

    a = 1 + random(y, 1, encoder.W-1)
    b = random(y, 2, encoder.W)

    if d < 4: 
        d1 = 2 + random(X, 3, 2)
    else: 
        d1 = 2

    a1 = 1 + random(X, 4, encoder.P1-1)
    b1 = random(X, 5, encoder.P1)

    return (d, a, b, d1, a1, b1)



def random(y, i, m):
    '''
    Pseudo random numbers as described here http://tools.ietf.org/html/rfc6330#section-5.3.5.1
    '''
    x0 = int( (y + i) % math.pow(2, 8) )
    x1 = int( (math.floor( y/math.pow(2, 8) ) + i) % math.pow(2, 8) )
    x2 = int( (math.floor( y/math.pow(2, 16) ) + i) % math.pow(2, 8) )
    x3 = int( (math.floor( y/math.pow(2, 24) ) + i) % math.pow(2, 8) )

    v0 = constants.V0[int(x0)]
    v1 = constants.V1[int(x1)]
    v2 = constants.V2[int(x2)]
    v3 = constants.V3[int(x3)]

    v01 = operator.xor(v0, v1)
    v012 = operator.xor(v01, v2)
    v0123 = operator.xor(v012, v3)
    ret_rand = v0123 % m
    return ret_rand

def degree(v, W):
    """
    Generate degree based on http://tools.ietf.org/html/rfc6330#section-5.3.5.2
    If v is out of range return None
    """
    if v >= 1048576:
        return None
    f = constants.f
    d = None
    for i in xrange(len(f)):
        if v < f[i]:
            d = i
            break

    if d < (W-2):
        return d
    return W-2


def ba_from_int(length, x):
    ba = bitarray(length)
    ba.setall(0)
    shift_x = x
    for i in xrange(len(ba)):
        ba[i] = shift_x % 2
        shift_x = shift_x >> 1

    return ba


def int_from_ba(ba):
    x = 0
    for i in xrange(len(ba)):
        x += int((ba[i] * math.pow(2, i)))

    return x
