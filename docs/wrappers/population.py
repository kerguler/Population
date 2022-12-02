import atexit
import numpy
from ctypes import *
import numpy.ctypeslib as npct
array_1d_double = npct.ndpointer(dtype=numpy.float64, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=numpy.int32, ndim=1, flags='CONTIGUOUS')

class model:
    def __init__(self, filename):
        self.filename = filename
        self.dylib = cdll.LoadLibrary(self.filename)
        #
        self.init = self.dylib.init
        self.init.restype = None
        self.init.argtypes = [array_1d_int, array_1d_int]
        self.numpar = numpy.arange(1,dtype=numpy.int32)
        self.nummet = numpy.arange(1,dtype=numpy.int32)
        ret = self.init(self.numpar,self.nummet)
        self.numpar = self.numpar[0]
        self.nummet = self.nummet[0]
        #
        atexit.register(self.dylib.destroy)
        #
        try:
            self.parnames = self.dylib.parnames
            self.parnames.restype = None
            self.parnames.argtypes = [POINTER(c_char_p)]
            temp = (c_char_p * (self.nummet+self.numpar))(256)
            ret = self.parnames(temp)
            temp = numpy.array([str(elm,'utf-8') for elm in temp])
            self.metnames = numpy.copy(temp[:self.nummet])
            self.parnames = numpy.copy(temp[-self.numpar:])
        except:
            print("Falling back to default parameters")
            self.metnames = numpy.array(["coln%d" %(n) for n in range(self.nummet)])
            self.parnames = numpy.array(["par%d" %(n) for n in range(self.numpar)])
        #
        self.metids = {}
        for elm in self.metnames:
            self.metids[elm] = numpy.where(elm==self.metnames)[0][0]
        #
        self.parids = {}
        for elm in self.parnames:
            self.parids[elm] = numpy.where(elm==self.parnames)[0][0]
        #
        self.csim = self.dylib.sim
        self.csim.restype = None
        self.csim.argtypes = [array_1d_double,
                              array_1d_double,
                              array_1d_int,
                              array_1d_int,
                              array_1d_double,
                              array_1d_int]
        #
    def sim(self,envir,pr,ftime,rep=1):
        """
            Note: Final time point is ftime - 1
        """
        envir = numpy.array(envir)
        pr = numpy.array(pr)
        ftime = numpy.array(ftime, dtype=numpy.int32, ndmin=1)
        rep = numpy.array(rep, dtype=numpy.int32, ndmin=1)
        ret = numpy.ndarray(rep*ftime*self.nummet, dtype=numpy.float64)
        success = numpy.array(0, dtype=numpy.int32, ndmin=1)
        self.csim(envir,
                  pr,
                  ftime,
                  rep,
                  ret,
                  success)
        ret = numpy.array(ret).reshape((rep[0],ftime[0],self.nummet))
        return ret

class parser:
    def __init__(self, json) -> None:
        self.json = json
        #
    def write_header(self):
        print("#include \"population.h\"")
        print()
        print("extern gsl_rng *RANDOM;")
        print()
        #
    def write_init(self):
        print("void init(int *np, int *nm) {")
        print("    spop2_random_init();")
        print()
        print("    *np = NumPar;")
        print("    *nm = NumMet;")
        print("}")
        print()
        #
    def write_destroy(self):
        print("void destroy(void) {")
        print("    spop2_random_destroy();")
        print("}")
        print()
        #
        
"""
TEST
"""
if __name__ == '__main__':
    pass