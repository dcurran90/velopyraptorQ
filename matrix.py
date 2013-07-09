import constants
from octet_array import octet_array

def divide_octet(u, v):
    if not u.count() or not v.count():
        return octet_array('00000000')
    u_int = u.val()
    v_int = v.val()
    log1 = constants.OCT_LOG[u_int]
    log2 = constants.OCT_LOG[v_int]
    val = constants.OCT_EXP[log1 - log2 + 255]

    return octet_array.from_val(val)

def multiply_octet(u, v):
    if not u.count() or not v.count():
        return octet_array('00000000')
    u_int = u.val()
    v_int = v.val()
    log1 = constants.OCT_LOG[u_int]
    log2 = constants.OCT_LOG[v_int]
    val = constants.OCT_EXP[log1 + log2]
    return octet_array.from_val(val)

def multiply_octet_by_row(beta, u):
    result = octet_array()
    for i in xrange(len(u)/8):
        val = multiply_octet(beta, u.get(i))
        result += val
    return result

def inverse_octet(u):
    if not u.count():
        return octet_array('00000000')
    u_int = u.val()
    val = constants.OCT_EXP[255 - constants.OCT_LOG[u_int]]
    return octet_array.from_val(val)

def octet_mat_multiply_bit(mat1, mat2):
    '''
    u*v col in u times rows in v
    '''
    if len(mat1[0])/8 != len(mat2):
        raise Exception("Cannot multiply matrices. Must be mxn and nxp to multiply.")

    new_mat = []
    for i in xrange(len(mat1)):
        ba = octet_array(len(mat2[0]))
        ba.setall(0)
        new_mat.append(ba)

    for i in xrange(len(mat1)):
        for j in xrange(len(mat2[0])/8):
            curr_val = octet_array('00000000')
            for k in xrange(len(mat2)):
                ba1 = mat1[i].get(k)
                ba2 = mat2[k].get(j)
                curr_val = curr_val ^ multiply_octet(ba1, ba2)
            new_mat[i].set(j, curr_val)

    return new_mat


    