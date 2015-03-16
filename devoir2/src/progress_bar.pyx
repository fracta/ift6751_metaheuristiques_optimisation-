"""simple progress bar"""

import sys
import time


cdef class ProgressBar:
    """simple progress bar class"""
    cdef str name
    cdef int bar_length

    def __init__(self, str name, int bar_length=25):
        self.name = name
        self.bar_length = bar_length

    cpdef update(self, float percentage):
        assert(0 <= percentage <= 1)
        block = int(round(self.bar_length*percentage))
        text = "\r[{0}] {1:.0f}% {2}".format( "#"*block + "-"*(self.bar_length-block), percentage*100, self.name)
        sys.stdout.write(text)
        sys.stdout.flush()

    cpdef clean(self):
        self.update(1.)
        sys.stdout.write("\n")
        sys.stdout.flush()

def test():
    # simple test
    print "progress : 0->1"
    p = ProgressBar("test")
    for i in range(100):
        time.sleep(0.1)
        p.update(i/100.0)
    return
