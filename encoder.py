from octet_array import octet_array
from raptor import RaptorQ

class Encoder(RaptorQ):
    def __init__(self, symbols):
        #Figure out sub-blocks. Basically split up "chunk" so it can fit in memory. Then split that into K pieces.
        #We'll just use one block for now
        super(Encoder, self).__init__(symbols)

        #Add padding symbols to get to KP
        for i in xrange(self.KP-self.K):
            ba = octet_array(symbols.symbolsize*8)
            ba.setall(0)
            symbols.append((self.K+i, ba))

        self.calculate_i_symbols()