from raptor import RaptorQ

class Decoder(RaptorQ):

    def __init__(self, symbols):
        #Figure out sub-blocks. Basically split up "chunk" so it can fit in memory. Then split that into K pieces.
        #We'll just use one block for now
        super(Decoder, self).__init__(symbols)

    def append(self, symbol):
        self.symbols.append(symbol)

    def decode(self):
        self.calculate_i_symbols()