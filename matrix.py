import copy
import math
import constants
import param_gen
from bitarray import bitarray

def divide_octet(u, v):
    if not u or not v:
        return 0
    log1 = constants.OCT_LOG[u]
    log2 = constants.OCT_LOG[v]
    val = constants.OCT_EXP[log1 - log2 + 255]
    return val

def multiply_octet(u, v):
    if not u or not v:
        return 0
    log1 = constants.OCT_LOG[u]
    log2 = constants.OCT_LOG[v]
    val = constants.OCT_EXP[log1 + log2]
    return val

def divide_octet_ba(u, v):
    if not u.count() or not v.count():
        return bitarray('00000000')
    u_int = param_gen.int_from_ba(u)
    v_int = param_gen.int_from_ba(v)
    log1 = constants.OCT_LOG[u_int]
    log2 = constants.OCT_LOG[v_int]
    val = constants.OCT_EXP[log1 - log2 + 255]
    return param_gen.ba_from_int(8, val)

def multiply_octet_ba(u, v):
    if not u.count() or not v.count():
        return bitarray('00000000')
    u_int = param_gen.int_from_ba(u)
    v_int = param_gen.int_from_ba(v)
    log1 = constants.OCT_LOG[u_int]
    log2 = constants.OCT_LOG[v_int]
    val = constants.OCT_EXP[log1 + log2]
    return param_gen.ba_from_int(8, val)

def multiply_octet_by_row(beta, u):
    result = bitarray()
    for i in xrange(len(u)/8):
        start = i*8
        end = (i+1)*8
        val = multiply_octet_ba(beta, u[start:end])
        result += val
    return result

def inverse_octet(u):
    if not u:
        return 0
    val = constants.OCT_EXP[255 - constants.OCT_LOG[u]]
    return val

def inverse_octet_ba(u):
    if not u.count():
        return bitarray('00000000')
    u_int = param_gen.int_from_ba(u)
    val = constants.OCT_EXP[255 - constants.OCT_LOG[u_int]]
    return param_gen.ba_from_int(8, val)

def xor_rows(u, v):
    new_row = []
    if len(u) != len(v):
        raise Exception("Cannot xor rows of different lengths")
    for i in xrange(len(u)):
        a = u[i]
        b = v[i]
        new_row.append(a^b)
    return new_row

def octet_mat_multiply_bit(mat1, mat2):
    '''
        u*v col in u times rows in v
    '''
    if len(mat1[0])/8 != len(mat2):
        raise Exception("Cannot multiply matrices. Must be mxn and nxp to multiply.")

    new_mat = []
    for i in xrange(len(mat1)):
        ba = bitarray(len(mat2[0]))
        ba.setall(0)
        new_mat.append(ba)

    for i in xrange(len(mat1)):
        for j in xrange(len(mat2[0])/8):
            start = j*8
            end = (j+1)*8
            curr_val = bitarray('00000000')
            for k in xrange(len(mat2)):
                ba1 = mat1[i][k*8:(k+1)*8]
                ba2 = mat2[k][start:end]
                curr_val = curr_val ^ multiply_octet_ba(ba1, ba2)
            new_mat[i][start:end] = curr_val

    return new_mat


def octet_mat_multiply(mat1, mat2):
    '''
    u*v col in u times rows in v
    '''

    if len(mat1[0]) != len(mat2):
        raise Exception("Cannot multiply matrices. Must be mxn and nxp to multiply.")

    new_mat = []
    for i in xrange(len(mat1)):
        ba_list = []
        for j in xrange(len(mat2[0])):
            ba = bitarray('00000000')
            ba_list.append(ba)
        new_mat.append(ba_list)

    for i in xrange(len(mat1)):
        for j in xrange(len(mat2[0])):
            curr_val = 0
            for k in xrange(len(mat2)):
                val1 = param_gen.int_from_ba(mat1[i][k])
                val2 = param_gen.int_from_ba(mat2[k][j])
                curr_val = curr_val ^ multiply_octet(val1, val2)

            new_mat[i][j] = param_gen.ba_from_int(8, curr_val)

    return new_mat


def octet_mat_multiply_int(mat1, mat2):
    '''
    u*v col in u times rows in v
    '''
    if len(mat1[0]) != len(mat2):
        raise Exception("Cannot multiply matrices. Must be mxn and nxp to multiply.")

    new_mat = []
    for i in xrange(len(mat1)):
        num_list = []
        for j in xrange(len(mat2[0])):
            num_list.append(0)
        new_mat.append(num_list)

    for i in xrange(len(mat1)):
        for j in xrange(len(mat2[0])):
            curr_val = 0
            for k in xrange(len(mat2)):
                val1 = mat1[i][k]
                val2 = mat2[k][j]
                curr_val = curr_val ^ multiply_octet(val1, val2)

            new_mat[i][j] = curr_val

    return new_mat

def identity_int(n):
    """
    Creates an n by n identity matrix
    Arguments:
    n -- Integer indicates rows and columns
    """

    matrix = []
    for i in xrange(n):
        row = [0] * n
        row[i] = 1
        matrix.append(row)
    return matrix

def inverse_int(mat):
    '''
    return the inverse of mat
    '''
    # Check dimensions.  Matrix must be square
    size = len(mat)
    if not (size == len(mat[0])):
        raise Exception("Tried to invert a %s by %s matrix. Matrix must be square" % (size, len(mat[0])))

    # We dont want to change a
    mat = copy.deepcopy(mat)

    id = identity_int(size)
    for i in xrange(size):
        mat[i] += id[i]

    # First produce the upper triangular matrix (1's in the diangonal
    # and 0's below
    for k in xrange(size):
        #Find max in current column
        column_vector = [x[k] for x in mat[k:]]
        i_max = column_vector.index(max(column_vector)) + k
        if mat[i_max][k] == 0:
            raise Exception("Cannot invert matrix because it is singular")

        if i_max != k:
            tmp_row = mat[k]
            mat[k] = mat[i_max]
            mat[i_max] = tmp_row

        for i in xrange(k, size):
            if mat[i][k] > 1:
                oct_inv = inverse_octet(mat[i][k])
                mat[i] = [multiply_octet(x, oct_inv) for x in mat[i]]
            if i > k and mat[i][k] == 1:
                mat[i] = xor_rows(mat[k], mat[i])

    # We should have 1s for the diagonal and all 0s
    for i in xrange(size):
        if mat[i][i] != 1:
            raise Exception("Tried to invert invertible matrix.  Matrix does not have ones in diagonal in upper triangular form")

    for k in xrange(size-1, -1, -1):
        if mat[k][k] != 1:
            oct_inv = inverse_octet(mat[k][k])
            mat[k] = [multiply_octet(x, oct_inv) for x in mat[k]]

        for i in xrange(k-1, -1, -1):
            if mat[i][k] > 1:
                oct_inv = inverse_octet(mat[i][k])
                mat[i] = [multiply_octet(x, oct_inv) for x in mat[i]]
            if i < k and mat[i][k] == 1:
                mat[i] = xor_rows(mat[k], mat[i])

    return [x[size:] for x in mat]