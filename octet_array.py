import math
from bitarray import bitarray

class octet_array(bitarray):
    def __init__(self, *args, **kwargs):
        super(octet_array, self).__init__(*args, **kwargs)

    @classmethod
    def from_val(cls, val):
        new_array = octet_array(8)
        new_array.setall(0)
        shift_val = val
        for i in xrange(len(new_array)):
            new_array[i] = shift_val % 2
            shift_val = shift_val >> 1
        return new_array
        
    def get(self, index):
        start = index*8
        end = (index+1)*8
        return octet_array(self[start:end])

    def set(self, index, val):
        if len(val) != 8:
            raise Exception("Value too long to add to octet_array")
        start = index*8
        end = (index+1)*8
        self[start:end] = val
        
    def val(self):
        my_val = 0
        for i in xrange(len(self)):
            my_val += int((self[i] * math.pow(2, i)))
        return my_val