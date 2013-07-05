import chunker
import encoder
import decoder
import raptor
from block import Source as EncodedBlock


K = 16 
symbolsize = 16

file_chunker = chunker.FileChunker(K, symbolsize, '../war_and_peace.txt')
block = file_chunker.chunk()

my_encoder = encoder.Encoder(K, block)

enc_symbols = my_encoder.enc_symbols
enc_block = EncodedBlock(K, symbolsize, 0)

for i in xrange(0, len(enc_symbols), 2):
#for i in xrange(2, my_encoder.KP+2):
    enc_block.append( (i, enc_symbols[i]) )

my_decoder = decoder.Decoder(enc_block)

decoded_sym = my_decoder.decoded_sym



for x in xrange(len(block)):
    if block[x] == decoded_sym[x]:
        print x, "OK"
    else:
        print x, "BAD"