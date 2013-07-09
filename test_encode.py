import chunker
import encoder
import decoder
import raptor
from block import Source as EncodedBlock


K = 16 
symbolsize = 16

file_chunker = chunker.FileChunker(K, symbolsize, '../war_and_peace.txt')
block = file_chunker.chunk()
for index in xrange(len(block)):
    val = block[index]
    block[index] = (index, val)




my_encoder = encoder.Encoder(block)

enc_symbols = [] #my_encoder.enc_symbols
for i in xrange(my_encoder.L):
    enc_symbols.append(my_encoder.next())

enc_block = EncodedBlock(K, symbolsize, 0)
for i in xrange(18, len(enc_symbols)):
#for i in xrange(0, my_encoder.KP):
    enc_block.append(enc_symbols[i])

my_decoder = decoder.Decoder(enc_block)
my_decoder.decode()


decoded_sym = []
for i in xrange(len(block)):
    sym_id, sym = my_decoder.next()
    decoded_sym.append(sym)

#decoded_sym = my_decoder.decoded_sym

for x in xrange(len(block)):
    if block[x][1] == decoded_sym[x]:
        print x, "OK"
    else:
        print x, "BAD"