from bitarray import bitarray

class octet_array(bitarray):
    def __init__(self, *args, **kwargs):
        super(octet_array, self).__init__(*args, **kwargs)

    def get_octet(self, index):
        start = index*8
        end = (index+1)*8
        return self[start:end]