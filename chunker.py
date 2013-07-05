from bitarray import bitarray
from block import Source as SourceBlock
import io

class FileChunker(object):

    def __init__(self, k, symbolsize, filename):

        """
        File chunker constructor

        Arguments:
        k -- Integer number of symbols per block
        symbol_size -- Integer Size of each symbol (IN BYTES)
        filename -- String name of file to chunk
        """

        self.block_id = 0
        self.k = k
        self.symbolsize = symbolsize # Bytes
        self.blocksize = self.symbolsize * self.k # Bytes
        
        try:
            self._f = io.open(filename, 'rb')
        except Exception:
            raise Exception('Unable to open file %s for reading.' % filename)

    def chunk(self):

        """
        Should return a block of bitarrays
        Attempts to read k symbols from the file and pads a block
        so that the block is k symbols of symbolsize * 8 bits
        """

        # Check to see file is still open
        if self._f.closed:
            return None
        
        block = SourceBlock(self.k, self.symbolsize, self.get_block_id())

        j = 0
        EOF = False

        while (j < self.k and not EOF):
            b = self._read()
            if b:
                block.append(b)
                j += 1
            else:
                EOF = True
                self.close()

        if len(block) == 0:
            return None

        block.pad()
        return block

    def _read(self):

        """
        Reads symbolsize bytes from the file
        Returns None if the length is 0
        Returns a bit array of symbolsize * 8 bits otherwise
        """

        # Read self.symbolsize bytes
        ba = bitarray()
        ba.frombytes(self._f.read(self.symbolsize))
        if len(ba) == 0:
            return None
        return ba

    def get_block_id(self):

        """
        Gets the current block id and increments to prepare for the next block.
        Returns the current block id
        """
        r = self.block_id
        self.block_id += 1
        return r

    def close(self):

        """
        Attempts to close the file associated with this chunker
        """

        try:
            self._f.close()
        except:
            pass
