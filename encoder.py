import constants
import param_gen
import matrix
import primes
import math
from bitarray import bitarray

class Encoder(object):
    def __init__(self, K, block):
        self.K = K
        self.block = block
        self.generate_params()

        #Add padding symbols to get to KP
        for i in xrange(self.KP-self.K):
            ba = bitarray(block.symbolsize*8)
            ba.setall(0)
            block.append(ba)

        #Figure out sub-blocks. Basically split up "chunk" so it can fit in memory. Then split that into K pieces.
        #We'll just use one block for now


        tuples = {
                    "d": [None] * self.L, 
                    "a": [None] * self.L, 
                    "b": [None] * self.L, 
                    "d1": [None] * self.L, 
                    "a1": [None] * self.L, 
                    "b1": [None] * self.L
                }

        for X in xrange(self.L):
            (tuples["d"][X], tuples["a"][X], tuples["b"][X], tuples["d1"][X], tuples["a1"][X], tuples["b1"][X]) = param_gen.tuples(self, X)

        #C[0], ..., C[B-1] denote the intermediate symbols that are LT symbols but not LDPC symbols.
        #C[B], ..., C[B+S-1] denote the S LDPC symbols that are also LT symbols.
        #C[W], ..., C[W+U-1] denote the intermediate symbols that are PI symbols but not HDPC symbols.
        #C[L-H], ..., C[L-1] denote the H HDPC symbols that are also PI symbols.

        #ALL MATRICES MUST BE OCTETS. ALL MATH ON THEM MUST BE GF(256)
        self.ldpc = self.gen_ldpc()
        self.hdpc, self.hdpc2 = self.gen_hdpc_IH()
        self.g_enc = self.gen_g_enc(tuples)
        

        #Get A from:
        #|ldpc    |
        #|hdpc2   | (this in int version of hdpc along with I_H)
        #|g_enc   |
        self.a = []
        for x in self.ldpc:
            row = []
            for y in x:
                row.append(int(y))
            self.a.append(row)

        for row in self.hdpc2:
            self.a.append(row)

        for x in self.g_enc:
            row = []
            for y in x:
                row.append(int(y))
            self.a.append(row)


        #Get inverse of A using terrible octet math. This is a LxL matrix of octets
        self.a_inverse = matrix.inverse_int(self.a)

        #D should be the symbols in a column. Split these into symbolsize pieces of octets
        #S+H padding symbols then K' source symbols
        self.d = []
        bitlength = self.block.symbolsize * 8
        for i in xrange(self.S + self.H):
            ba = bitarray(bitlength)
            ba.setall(0)
            self.d.append(ba)

        for symbol in block:
            self.d.append(symbol)

        for i in xrange(len(self.d)):
            symbol = self.d[i]
            self.d[i] = []
            for j in xrange(len(symbol)/8):
                start = 8*j
                end = start+8
                num = param_gen.int_from_ba(bitarray(symbol[start:end]))
                self.d[i].append(num)

        #multiply a_inverse and d to get the intermediate symbols (hopefully)
        #c_int will be the intermediate symbols as integers.
        #C is intermediates as bitarrays of length (symbolsize*8)
        self.c_int = matrix.octet_mat_multiply_int(self.a_inverse, self.d)

        self.C = []
        for x in self.c_int:
            row = bitarray()
            for y in x:
                row = row+param_gen.ba_from_int(8,y)
            self.C.append(row)

        self.enc_symbols = []
        for i in xrange(self.L):
            self.enc_symbols.append(self.encode(i, tuples))




    def encode(self, X, tuples):
        d = tuples["d"][X]
        a = tuples["a"][X]
        b = tuples["b"][X]
        d1 = tuples["d1"][X]
        a1 = tuples["a1"][X]
        b1 = tuples["b1"][X]

        result = self.C[b]
        for j in xrange(1,d):
            b = (b + a) % self.W
            result = result ^ self.C[b]

        while b1 >= self.P:
            b1 = (b1+a1) % self.P1

        result = result ^ self.C[self.W+b1]
        for j in xrange(1, d1):
            b1 = (b1 + a1) % self.P1
            while b1 >= self.P:
                b1 = (b1+a1) % self.P1
            result = result ^ self.C[self.W+b1]

        return result        


    def gen_g_enc(self, tuples):
        g_enc = []
        for i in xrange(self.KP):
            ba = bitarray(self.L)
            ba.setall(0)            

            d = tuples["d"][i]
            a = tuples["a"][i]
            b = tuples["b"][i]
            d1 = tuples["d1"][i]
            a1 = tuples["a1"][i]
            b1 = tuples["b1"][i]

            ba[b] = True
            for j in xrange(1,d):
                b = (b + a) % self.W
                ba[b] = True

            while b1 >= self.P:
                b1 = (b1+a1) % self.P1

            ba[self.W + b1] = True
            for j in xrange(1, d1):
                b1 = (b1 + a1) % self.P1
                while b1 >= self.P:
                    b1 = (b1+a1) % self.P1
                ba[self.W+b1] = True

            g_enc.append(ba)
        return g_enc


    def gen_ldpc(self):
        ldpc = []
        for i in xrange(self.S):
            ba = bitarray(self.L)
            ba.setall(0)
            ldpc.append(ba)

        for i in xrange(self.B):
            a = int(1 + math.floor(i / self.S))
            b = i % self.S
            ldpc[b][i] = True
            b = (b+a) % self.S
            ldpc[b][i] = True
            b = (b+a) % self.S
            ldpc[b][i] = True

        for i in xrange(self.S):
            row = i
            col = i + self.B
            ldpc[row][col] = True

        for i in xrange(self.S):
            a = i % self.P
            b = (i+1) % self.P
            ldpc[i][self.W + a] = True
            ldpc[i][self.W + b] = True

        return ldpc


    def gen_hdpc_IH(self):
        MT = []
        MT2 = []
        for i in xrange(self.H):
            ba_list = []
            num_list = []
            for j in xrange(self.KP + self.S):
                ba = bitarray('00000000')
                ba_list.append(ba)
                num_list.append(0)
            MT.append(ba_list)
            MT2.append(num_list)


        for i in xrange(self.H):
            for j in xrange(self.KP + self.S):
                rand_1 = param_gen.random(j+1, 6, self.H)
                rand_2 = (param_gen.random(j+1,6,self.H) + param_gen.random(j+1,7,self.H-1) + 1) % self.H
                if i == rand_1 or i == rand_2:
                    MT[i][j] = bitarray('10000000')
                    MT2[i][j] = 1

        j = self.KP + self.S - 1
        for i in xrange(self.H):
            alpha_exp = param_gen.ba_from_int(8, constants.OCT_EXP[i])
            MT[i][j] = alpha_exp
            MT2[i][j] = constants.OCT_EXP[i]

        
        GAMMA = []
        GAMMA2 = []
        for i in xrange(self.KP + self.S):
            ba_list = []
            num_list = []
            for j in xrange(self.KP + self.S):
                ba = bitarray('00000000')
                num = 0
                if i >= j:
                    ba = param_gen.ba_from_int(8, constants.OCT_EXP[i-j])
                    num = constants.OCT_EXP[i-j]
                ba_list.append(ba)
                num_list.append(num)
            GAMMA.append(ba_list)
            GAMMA2.append(num_list)

        hdpc = matrix.octet_mat_multiply(MT, GAMMA)
        hdpc2 = matrix.octet_mat_multiply_int(MT2, GAMMA2)

        for i in xrange(self.H):
            for j in xrange(self.H):
                if i == j:
                    ba = param_gen.ba_from_int(8, 1)
                    num = 1
                else:
                    ba = bitarray('00000000')
                    num = 0
                hdpc[i].append(ba)
                hdpc2[i].append(num)

        return hdpc, hdpc2



    def generate_params(self):
        #C[0]-C[B] are intermediate
        #C[B]-C[B+S-1] are LDPC
        #C[W]-C[W+U-1] intermediate PI symbols that aren't HDPC
        #C[L-H]-C[L-1] intermediate PI symbols that are HDPC
        #Tuples are (K',J(K'),S(K'),H(K'),W(K'))
        for (self.KP, self.J, self.S, self.H, self.W) in constants.sys_indices:
            if self.KP >= self.K:
                break
                
        self.L = self.KP + self.S + self.H
        self.P = self.L - self.W
        self.P1 = primes.next(self.P)
        self.U = self.P - self.H
        self.B = self.W - self.S
