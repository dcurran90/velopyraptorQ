import constants
import param_gen
import networkx
import matrix
import primes
import math
import copy
from octet_array import octet_array
from schedule import Schedule

class RaptorQ(object):
    def __init__(self, symbols):
        self.K = symbols.k
        self.symbols = symbols
        #Tuples are (K',J(K'),S(K'),H(K'),W(K'))        
        for (self.KP, self.J, self.S, self.H, self.W) in constants.sys_indices:
            if self.KP >= self.K:
                break

        self.L = self.KP + self.S + self.H
        self.P = self.L - self.W
        self.P1 = primes.next(self.P)
        self.U = self.P - self.H
        self.B = self.W - self.S

        self.current_id = 0

    def next(self):
        sym_id, ret_sym = self.current_id, self.encode(self.current_id)
        self.current_id += 1
        return sym_id, ret_sym

    def encode(self, X):
        d, a, b, d1, a1, b1 = param_gen.tuples(self, X)

        result = copy.deepcopy(self.i_symbols[b])
        for j in xrange(1,d):
            b = (b + a) % self.W
            result = result ^ self.i_symbols[b]

        while b1 >= self.P:
            b1 = (b1+a1) % self.P1
        result = result ^ self.i_symbols[self.W+b1]
        
        for j in xrange(1, d1):
            b1 = (b1 + a1) % self.P1
            while b1 >= self.P:
                b1 = (b1+a1) % self.P1
            result = result ^ self.i_symbols[self.W+b1]

        return result

    def calculate_i_symbols(self):
        self.a = self.gen_a()

        schedule = self.create_schedule()

        D = self.calculate_d()
        
        for beta, xor_row, target_row in schedule.ops:
            new_row = matrix.multiply_octet_by_row(beta, D[xor_row])
            if target_row != None:
                D[target_row] ^= new_row
            else:
                D[xor_row] = new_row

        self.i_symbols = [None]*self.L
        for i in xrange(self.L):
            self.i_symbols[schedule.c[i]] = D[schedule.d[i]]

    def calculate_d(self):
        d = []
        symbolsize = self.symbols.symbolsize*8

        for i in xrange(self.S+self.H):
            ba = octet_array(symbolsize)
            ba.setall(0)
            d.append(ba)

        for i, symbol in self.symbols:
            d.append(octet_array(copy.deepcopy(symbol)))

        return d

    def create_schedule(self):
        n = len(self.symbols)
        m = self.S + self.H + n
        self.X = copy.deepcopy(self.a)

        i = 0
        u = self.P
        c = range(self.L)
        d = range(m)
        self.X = copy.deepcopy(self.a)


        schedule = Schedule(self.L, m)
        schedule_X = Schedule(self.L, m)


        #CHECK OUT RIGHT ABOVE SECTION 5.4.2.3
        #EFFICIENCY CAN IMPROVE WHEN YOU DONT APPLY ROW OPS TILL THEY"RE CHOSEN IN DECODER
        while (i + u) < self.L:
            r, rows_with_r = self.rows_in_v_with_min_r(self.a, m, i, u)
            if r == 0:
                raise Exception("Unable to decode.  No nonzero row to choose from v")

            #Separate hdpc and non_hdpc rows. Non hdpc should always be handled first
            non_hdpc_rows = []
            hdpc_rows = []
            for row_index in rows_with_r:
                if row_index<self.S or row_index>=(self.S+self.H):
                    non_hdpc_rows.append(row_index)
                else:
                    hdpc_rows.append(row_index)


            # If r is two and we have non_hdpc rows then use the graph to find row, otherwise choose any row
            # If r is not two choose row with minimum original degree
            if r == 2:
                if len(non_hdpc_rows) > 0:
                    chosen_row_index = self.choose_row_from_graph(self.a, m, i, u, non_hdpc_rows)
                else:
                    chosen_row_index = rows_with_r[0]
            else:
                if len(non_hdpc_rows) > 0:
                    chosen_row_index = self.choose_min_degree_row(self.a, m, i, u, non_hdpc_rows)
                else:
                    chosen_row_index = self.choose_min_degree_row(self.a, m, i, u, rows_with_r)

            self.exchange_row(self.a, i, chosen_row_index, schedule)
            self.exchange_row(self.X, i, chosen_row_index, schedule_X)

            #Reorder columns, 1 non-zero followed by 0's then r-1 non-zeros
            ba = self.a[i].get(i)
            if ba.count() == 0:
                for col in xrange(self.L-u):
                    ba = self.a[i].get(col)
                    if ba.count() > 0:
                        self.exchange_column(self.a, i, col, schedule)
                        self.exchange_column(self.X, i, col, schedule_X)

            back_col = (self.L-u)-1
            front_col = i+1
            while front_col < back_col:
                front_ba = self.a[i].get(front_col)
                while front_ba.count() == 0 and front_col < back_col:
                    front_col += 1
                    front_ba = self.a[i].get(front_col)

                back_ba = self.a[i].get(back_col)
                while back_ba.count() > 0 and back_col > front_col:
                    back_col -= 1
                    back_ba = self.a[i].get(back_col)

                if front_col < back_col:
                    self.exchange_column(self.a, front_col, back_col, schedule)
                    self.exchange_column(self.X, front_col, back_col, schedule_X)

            chosen_entry = self.a[i].get(i)
            chosen_val = chosen_entry.val()#param_gen.int_from_ba(chosen_entry)
            for row_index in xrange(i+1, m):
                curr_entry = self.a[row_index].get(i)
                curr_val = curr_entry.val()#param_gen.int_from_ba(curr_entry)
                if curr_val != 0:
                    beta = matrix.divide_octet(curr_entry, chosen_entry)
                    self.xor_row(self.a, beta, i, row_index, schedule)

            if chosen_val > 1:
                beta = matrix.divide_octet(octet_array('10000000'), chosen_entry)
                self.multiply_row(self.a, beta, i, schedule)

            i += 1
            u += (r-1)

        #Discard rows/columns of X
        #So we're i x i matrix in lower triangular form
        del self.X[i:]
        for row in self.X:
            del row[i*8:]

        #GAUSSIAN ON U_LOWER
        for col_index in xrange(i, self.L):
            pivot_ba = self.a[col_index].get(col_index)
            pivot_val = pivot_ba.val()# param_gen.int_from_ba(pivot_ba)
            if pivot_val == 0:
                for row_index in xrange(col_index+1, m):
                    pivot_ba = self.a[row_index].get(col_index)
                    pivot_val = pivot_ba.val()#param_gen.int_from_ba(pivot_ba)
                    if pivot_val > 0:
                        self.exchange_row(self.a, col_index, row_index, schedule)
                if pivot_val == 0:
                    raise Exception("U lower is of less rank than %s." % u)

            if pivot_val != 1:
                entry_inv = matrix.inverse_octet(pivot_ba)
                self.multiply_row(self.a, entry_inv, col_index, schedule)
                pivot_ba = self.a[col_index].get(col_index)

            for row_index in xrange(col_index+1, m):
                curr_entry = self.a[row_index][col_index*8: (col_index+1)*8]
                curr_val = curr_entry.val()#param_gen.int_from_ba(curr_entry)
                if curr_val != 0:
                    beta = octet_array('10000000')
                    if curr_val > 1:
                        beta = matrix.divide_octet(curr_entry, pivot_ba)
                    self.xor_row(self.a, beta, col_index, row_index, schedule)

        # U Lower should now be in upper triangular form. now attack the top
        for col_index in xrange(self.L - 1, self.L - u - 1, -1):
            pivot_ba = self.a[col_index].get(col_index)
            pivot_val = pivot_ba.val()#param_gen.int_from_ba(pivot_ba)
            for row_index in xrange(i, col_index):
                curr_ba = self.a[row_index].get(col_index)
                curr_val = curr_ba.val()#param_gen.int_from_ba(curr_ba)
                if curr_val != 0:
                    beta = octet_array('10000000')
                    if curr_val > 1:
                        beta = matrix.divide_octet(curr_ba, pivot_ba)
                    self.xor_row(self.a, beta, col_index, row_index, schedule)
        #GAUSSIAN ON U_LOWER DONE

        # Discard any rows left after l
        del self.a[self.L:]

        #Multiply X by A to make U_Upper sparse. this means top-left of A is X
        #new_top_a = matrix.octet_mat_multiply_bit(self.X, self.a[:i])
        #for row_index in xrange(len(new_top_a)):
        #    self.a[row_index] = new_top_a[row_index]

        #Zero out U_Upper with multiples of the identity
        for row_index in xrange(i):
            for col_index in xrange(i, self.L):
                curr_ba = self.a[row_index].get(col_index)
                curr_val = curr_ba.val()#param_gen.int_from_ba(curr_ba)
                if curr_val != 0:
                    self.xor_row(self.a, curr_ba, col_index, row_index, schedule)

        #Make the rest of A be the identity matrix
        '''
        for j in xrange(0, i):
            curr_ba = self.a[j][j*8:(j+1)*8]
            curr_val = curr_ba.val()#param_gen.int_from_ba(curr_ba)
            if curr_val != 1:
                curr_inv = matrix.inverse_octet_ba(curr_ba)
                self.multiply_row(self.a, curr_inv, j, schedule)
                curr_ba = self.a[j][j*8:(j+1)*8]
            for l in xrange(0, j):
                curr_ba = self.a[j][l*8:(l+1)*8]
                curr_val = curr_ba.val()#param_gen.int_from_ba(curr_ba)
                if curr_val != 0:
                    pivot_ba = self.a[l][l*8:(l+1)*8]
                    self.xor_row(self.a, curr_ba, l, j, schedule)
        '''

        return schedule

    @classmethod
    def xor_row(cls, a, beta, r1, r2, schedule):

        """
        XORS r2 into r1 and records the operation
        within the schedule

        Arguments:
        a -- List of l sized bitarrays
        beta -- bitarray indicating multiplier with row r1
        r1 -- Integer source row id
        r2 -- Integer target row id
        schedule -- Schedule to record the operation in
        """

        
        # XOR r2 of a into r1 of a
        new_row = matrix.multiply_octet_by_row(beta, a[r1])
        a[r2] ^= new_row

        # Schedule the xor
        schedule.xor(beta, r1, r2)

    @classmethod
    def multiply_row(cls, a, beta, r1, schedule):

        """
        Multiplies A[r1] by beta and records in schedule
        Beta must be bitarray

        Arguments:
        beta -- Integer indicating multiplier with row r1
        r1 -- Integer indicating target row
        """
        a[r1] = matrix.multiply_octet_by_row(beta, a[r1])

        #Schedule the multiply
        schedule.multiply(beta, r1)


    @classmethod
    def exchange_column(cls, a, c1, c2, schedule):

        """
        Exchanges column c1 of a with column c2 of a and records the operation
        in the schedule

        Arguments:
        a -- List of bit arrays representing a
        c1 -- Integer first column id
        c2 -- Integer second column id
        schedule -- Schedule of operations performed upon a        
        """

        # Exchange the columns c1 and c2 in a
        for i in xrange(len(a)):
            temp = a[i].get(c1)
            a[i].set(c1, a[i].get(c2))
            a[i].set(c2, temp)

        # Record the operation
        schedule.exchange_column(c1, c2)

    @classmethod
    def exchange_row(cls, a, r1, r2, schedule):

        """
        Exchanges row r1 of a with row r2 of a and records the operation
        in the schedule

        Arguments:
        a -- list of bitarrays representing a
        r1 -- Integer id of first row to exchange
        r2 -- Integer id of second row to exchange
        schedule -- Schedule to record the operation in
        """

        # Exchange r1 with r2 of a
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp

        # Record the operation
        schedule.exchange_row(r1, r2)


    def choose_min_degree_row(self, a, m, i, u, rows_with_r):

        """
        Chooses a minimum degree row out of rows with r

        Arguments:
        a -- List of bitarrays representing matrix a
        m -- Integer n + s + h(a should have m rows)
        i -- Integer representing the i'th iteration in reducing matrix V
        u -- Integer representing number of columns in matix U
        rows_with_r -- List of row columns sharing the same number of ones
            in matrix V
        """
        # Calculate number of non-zero entries in columns of amongst V.
        all_degrees = {}
        for column in xrange(i, self.L - u):
            all_degrees[column] = []
            for row in xrange(i, m):
                octet_ba = a[row].get(column)
                if octet_ba.count() > 0:
                    all_degrees[column].append(row)

        # Degrees contains only columns that contain one of the rows in question
        degrees = {}
        for row in rows_with_r:
            for col in all_degrees:
                if row in all_degrees[col]:
                    degrees[col] = all_degrees[col]

        # Find minimum column
        min_degree = m + 1
        min_column = None
        for column in degrees:
            if not (len(degrees[column]) == 0):
                if len(degrees[column]) < min_degree:
                    min_degree = len(degrees[column])
                    min_column = degrees[column]


        # Return first row of min column that is among rows_with_r
        for row in min_column:
            if row in rows_with_r:
                return row


    def choose_row_from_graph(self, a, m, i, u, graph_rows):

        """
        Builds a graph from rows where rows are edges and columns are vertices.
        Then chooses the first edge from the largest component

        Arguments:
        a -- List of bitarrays representing matrix A
        m -- Integer total number of rows in A
        i -- Integer representing i'th iteration of reducing V
        u -- Integer representing number of columns in u
        rows_with_r -- List of row indexes that share the same number of ones
            in V
        """
        #find rows with exactly 2 non-zero entries
        graph = networkx.Graph()
        for row in graph_rows:
            vertices = []
            for vertex in xrange(i, self.L - u):
                if a[row][vertex*8]:
                    vertices.append(vertex)
            v1, v2 = tuple(vertices)
            graph.add_edge(v1, v2, row_index=row)

        # Calculate components in graph
        components = networkx.connected_component_subgraphs(graph)

        # Find the max component
        max_component = None
        max_size = 0
        for c in components:
            edges = c.edges(data=True)
            if len(edges) > max_size:
                max_component = edges
                max_size = len(edges)

        _, _, data = max_component[0]
        row = data['row_index']
        return row


    def rows_in_v_with_min_r(self, a, m, i, u):

        """
        Returns a tuple with the minimum number of 1s in a row in v
        and the indexes of the rows containing that number of 1s

        Arguments:
        a -- List of bitarrays representing the matrix A
        m -- Integer total number of rows in A
        i -- Integer indicating i'th iteration of reducing V in A
        u -- Integer number of columns in matrix U

        Returns tuple (minimum r, list of rows with minimum r)
        """

        # Find minimum number of ones in a row in sub matrix v.
        min_r = None
        rows_with_min_r = []
        for row_index in xrange(i, m):
            v_row = a[row_index]

            # let r be the number of ones in a row in v
            r = 0
            for col_index in xrange(i, self.L-u):
                ba = v_row.get(col_index)
                if ba.count() > 0:
                    r += 1

            # Ignore rows 
            if r == 0:
                continue

            if min_r is None or r < min_r:
                min_r = r
                rows_with_min_r = [row_index]
            elif r == min_r:
                rows_with_min_r.append(row_index)

        return min_r, rows_with_min_r


    

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
        KP  |g_enc                |  number of symbols we're given)
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
        for i, _ in self.symbols:
            ba = octet_array(self.L*8)
            ba.setall(0)

            d, a, b, d1, a1, b1 = param_gen.tuples(self, i)

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
            ba = octet_array()
            for j in xrange(self.KP + self.S):
                ba += octet_array('00000000')
            MT.append(ba)

        #If is is rand_1 or rand_2 set values to 1
        for i in xrange(self.H):
            for j in xrange(self.KP + self.S):
                start = j*8
                rand_1 = param_gen.random(j+1,6,self.H)
                rand_2 = (rand_1 + param_gen.random(j+1,7,self.H-1) + 1) % self.H
                if i == rand_1 or i == rand_2:
                    MT[i][start] = True

        #For column K'+S-1 set all row values to alpha^^row_index
        j = self.KP + self.S - 1
        start = j*8
        end = (j+1)*8
        for i in xrange(self.H):
            alpha_exp = octet_array.from_val(constants.OCT_EXP[i])
            MT[i][start:end] = alpha_exp

        #GAMMA[i,j] = alpha ^^ (i-j) for i >= j, 0 otherwise
        GAMMA = []
        for i in xrange(self.KP + self.S):
            ba = octet_array()
            for j in xrange(self.KP + self.S):
                if i >= j:
                    ba += octet_array.from_val(constants.OCT_EXP[i-j])
                else:
                    ba += octet_array('00000000')
            GAMMA.append(ba)

        hdpc = matrix.octet_mat_multiply_bit(MT, GAMMA)

        #Throw I_H on the end
        for i in xrange(self.H):
            for j in xrange(self.H):
                if i == j:
                    hdpc[i] += octet_array('10000000')
                else:
                    hdpc[i] += octet_array('00000000')

        return hdpc

    def gen_ldpc(self):
        '''
        LDPC is three parts really
        LDPC,1 and 2 have 1's in places correspoding to indices of C
        that are added based on LDPC relations
        '''
        ldpc = []
        for i in xrange(self.S):
            ba = octet_array(self.L*8)
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
