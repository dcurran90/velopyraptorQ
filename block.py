from bitarray import bitarray

class Source(list):

    def __init__(self, k, symbolsize, block_id):

        """
        Block constructor
        
        Arguments:
        k -- number of symbols
        symbolsize -- Size of each symbol in bytes
        block_id -- id of this block
        """
        self.k = k
        self.symbolsize = symbolsize
        self.id = block_id
        self.padding = 0 # Bytes

    def pad(self):

        """
        Pads a block to have k bitarrays of symbolsize * 8 bits length
        
        Arguments:
        k -- number of bit arrays
        symbolsize -- Size of each array in bytes
        """

        bitlength = self.symbolsize * 8

        for row in self:
            if len(row) < bitlength:
                self.padding += (bitlength - len(row)) / 8
                extension = bitarray(bitlength - len(row)) 
                extension.setall(False)
                row.extend(extension)
                # Last provided row should be the only row that needs row padding
                break

        if len(self) < self.k:
            for i in xrange(len(self), self.k):
                padrow = bitarray(bitlength)
                padrow.setall(False)
                self.append(padrow)
                self.padding += bitlength / 8
