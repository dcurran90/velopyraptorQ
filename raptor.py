import constants
import param_gen
import matrix
import primes
import math
from bitarray import bitarray

class RaptorQ(object):
    def __init__(self, block):
        self.K = block.k
        self.block = block
        self.generate_params()
        self.a = self.gen_a()

    def gen_a(self):
        '''
        Produces matrix A which is used to create intermediate symbols.
        A look like the following
            +---------------------+
        S   |ldpc,1 | I_S | ldpc,2| 
            +---------------------+
        H   |hdpc           | I_H |
            +---------------------+
            |                     | (Should have height equal to 
        KP  |g_enc                |  number of blocks we're given)
            +---------------------+
        '''
        #These have the identity matries built into 
        #them so all we have to do is put them together
        ldpc = self.gen_ldpc()
        hdpc = self.gen_hdpc_IH()
        g_enc = self.gen_g_enc()
        return (ldpc + hdpc + g_enc)


    def gen_g_enc(self):
        '''
        Creates the g_enc matrix. Bottom part of A.
        Basically, use the tuples to find indices of C values that 
        will be xor'd together during encoding for each value of i (ESI).
        i corresponds to rows and the indices correspond to columns
        '''
        g_enc = []
        for i, _ in self.block:
            ba = bitarray(self.L*8)
            ba.setall(0)

            d, a, b, d1, a1, b1 = param_gen.tuples(self, i)
            #d = self.tuples["d"][i]
            #a = self.tuples["a"][i]
            #b = self.tuples["b"][i]
            #d1 = self.tuples["d1"][i]
            #a1 = self.tuples["a1"][i]
            #b1 = self.tuples["b1"][i]

            ba[b*8] = True
            for j in xrange(1,d):
                b = (b + a) % self.W
                ba[b*8] = True

            while b1 >= self.P:
                b1 = (b1+a1) % self.P1

            ba[(self.W + b1)*8] = True
            for j in xrange(1, d1):
                b1 = (b1 + a1) % self.P1
                while b1 >= self.P:
                    b1 = (b1+a1) % self.P1
                ba[(self.W+b1)*8] = True

            g_enc.append(ba)

        return g_enc

    def gen_hdpc_IH(self):
        '''
        HDPC is MT*GAMMA. 

        MT is Hx(K'+S) matrix made up of 1's for random values generated based on our columns
        or alpha^^i for column K'+S-1 for all rows

        GAMMA is (K'+S)x(K'+S) matrix of alpha^^i-j for i>=j or 0

        Append I_H to the end
        '''
        MT = []
        #Fill up MT with 0's
        for i in xrange(self.H):
            ba = bitarray()
            for j in xrange(self.KP + self.S):
                ba += bitarray('00000000')
            MT.append(ba)

        #If is is rand_1 or rand_2 set values to 1
        for i in xrange(self.H):
            for j in xrange(self.KP + self.S):
                start = j*8
                rand_1 = param_gen.random(j+1,6,self.H)
                rand_2 = (rand_1 + param_gen.random(j+1,7,self.H-1) + 1) % self.H
                if i == rand_1 or i == rand_2:
                    MT[i][start] = bitarray('1')

        #For column K'+S-1 set all row values to alpha^^row_index
        j = self.KP + self.S - 1
        start = j*8
        end = (j+1)*8
        for i in xrange(self.H):
            alpha_exp = param_gen.ba_from_int(8, constants.OCT_EXP[i])
            MT[i][start:end] = alpha_exp

        #GAMMA[i,j] = alpha ^^ (i-j) for i >= j, 0 otherwise
        GAMMA = []
        for i in xrange(self.KP + self.S):
            ba = bitarray()
            for j in xrange(self.KP + self.S):
                if i >= j:
                    ba += param_gen.ba_from_int(8, constants.OCT_EXP[i-j])
                else:
                    ba += bitarray('00000000')
            GAMMA.append(ba)

        hdpc = matrix.octet_mat_multiply_bit(MT, GAMMA)

        #Throw I_H on the end
        for i in xrange(self.H):
            for j in xrange(self.H):
                if i == j:
                    hdpc[i] += bitarray('10000000')
                else:
                    hdpc[i] += bitarray('00000000')

        return hdpc

    def gen_ldpc(self):
        '''
        LDPC is three parts really
        LDPC,1 and 2 have 1's in places correspoding to indices of C
        that are added based on LDPC relations
        '''
        ldpc = []
        for i in xrange(self.S):
            ba = bitarray(self.L*8)
            ba.setall(0)
            ldpc.append(ba)

        #Create LDPC,1
        for i in xrange(self.B):
            a = int(1 + math.floor(i / self.S))
            b = i % self.S
            ldpc[b][i*8] = True
            b = (b+a) % self.S
            ldpc[b][i*8] = True
            b = (b+a) % self.S
            ldpc[b][i*8] = True

        #Throw I_S in the middle of LDPC
        for i in xrange(self.S):
            row = i
            col = i + self.B
            ldpc[row][col*8] = True

        #Create LDPC,2
        for i in xrange(self.S):
            a = i % self.P
            b = (i+1) % self.P
            ldpc[i][(self.W + a)*8] = True
            ldpc[i][(self.W + b)*8] = True

        return ldpc

    def generate_params(self):
        '''
        Generate a bunch of parameters that are common to most functions
        TODO: ADD THE DESCRIPTIONS OF PARAMETERS
        '''
        #Tuples are (K',J(K'),S(K'),H(K'),W(K'))
        for (self.KP, self.J, self.S, self.H, self.W) in constants.sys_indices:
            if self.KP >= self.K:
                break
                
        self.L = self.KP + self.S + self.H
        self.P = self.L - self.W
        self.P1 = primes.next(self.P)
        self.U = self.P - self.H
        self.B = self.W - self.S

        self.tuples = {
                        "d": [None] * len(self.block),
                        "a": [None] * len(self.block),
                        "b": [None] * len(self.block),
                        "d1": [None] * len(self.block),
                        "a1": [None] * len(self.block),
                        "b1": [None] * len(self.block)
                    }
        for X in xrange(len(self.block)):
            (self.tuples["d"][X], self.tuples["a"][X], self.tuples["b"][X], self.tuples["d1"][X], self.tuples["a1"][X], self.tuples["b1"][X]) = param_gen.tuples(self, X)
