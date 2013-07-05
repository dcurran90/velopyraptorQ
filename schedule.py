class Schedule(object):

    """
    Class is used to keep track of operations on a matrix that can
    be mirrored onto another matrix
    """

    def __init__(self, l, m):

        """
        Initializes the sequences c and d

        Arguments:
        l -- Integer combined h + s + k(number of source symbols)
        m -- Integer combined h + s + n(number of encoding symbols.  n should be >= k)
        """

        # Let c[0] = 0, c[1] = 1,...,c[L-1] = L-1
        self.c = range(l)

        # Let d[0] = 0, d[1] = 1,...,d[M-1] = M-1
        self.d = range(m)

        # Init ops to empty list
        self.ops = []

    def xor(self, beta, r1, r2):

        """
        Indicates beta*A[r1] is xored into A[r2].  Appends a tuple(beta, d[r1], d[r2])
        to the list of ops
        Ordering of r1, r2 DOES matter

        Arguments:
        beta -- Integer indicating multiplier with row r1
        r1 -- Integer indicating source row
        r2 -- Integer indicating target row
        """
        self.ops.append((beta, self.d[r1], self.d[r2]))

    def multiply(self, beta, r1):

        """
        Indicates A[r1] is multiplied by beta.  Appends a tuple(beta, d[r1]) to the list of ops

        Arguments:
        beta -- Integer indicating multiplier with row r1
        r1 -- Integer indicating target row
        r2 -- Integer indicating source row
        """
        self.ops.append((beta, self.d[r1], None))        

    def exchange_row(self, r1, r2):

        """
        Indicates row r1 is swapped with row r2
        Ordering of r1, r2 DOES NOT matter

        Arguments:
        r1 -- Integer first row to swap
        r2 -- Integer second row to swap
        """
        
        swap = self.d[r1]
        self.d[r1] = self.d[r2]
        self.d[r2] = swap

    def exchange_column(self, c1, c2):

        """
        Indicates column c1 is swapped with column c2.
        Ordering of c1, c2 DOES NOT matter.

        Arguments:
        c1 -- Integer first column to swap
        c2 -- Integer second column to swap
        """
        swap = self.c[c1]
        self.c[c1] = self.c[c2]
        self.c[c2] = swap
