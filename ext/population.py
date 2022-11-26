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
        self.init.argtypes = [array_1d_int, array_1d_int, array_1d_int]
        self.numenv = numpy.arange(1,dtype=numpy.int32)
        self.numpar = numpy.arange(1,dtype=numpy.int32)
        self.nummet = numpy.arange(1,dtype=numpy.int32)
        ret = self.init(self.numenv,self.numpar,self.nummet)
        self.numenv = self.numenv[0]
        self.numpar = self.numpar[0]
        self.nummet = self.nummet[0]
        #
        atexit.register(self.dylib.destroy)
        #
        try:
            self.parameters = self.dylib.parameters
            self.parameters.restype = None
            self.parameters.argtypes = [POINTER(c_char_p)]
            temp = (c_char_p * (self.nummet+self.numpar))(256)
            ret = self.parameters(temp)
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
        csim = model.sim
        csim.restype = None
        csim.argtypes = [array_1d_double,
                         array_1d_double,
                         array_1d_int,
                         array_1d_double,
                         array_1d_int]
        #
    def sim(self,envir,pr):
        pr = numpy.array(pr)
        tf = temp.shape[0] + 1
        ret = numpy.ndarray(tf*8,dtype=numpy.float64)
        csim(tf,envirflat,pr,init,thr,ret)
        ret = ret.reshape((tf,8))
        return ret
    #
    def sim(self,clim,pr,cpr=[]):
        """
        Main simulation routine
        """
        fT = numpy.array(len(clim['dates']),dtype=numpy.int32,ndmin=1)
        envar = numpy.array(clim['envar'],dtype=numpy.float64)
        param = numpy.array(pr,dtype=numpy.float64)
        param = numpy.hstack([param,cpr])
        control = numpy.array(len(cpr)!=0,dtype=numpy.int32,ndmin=1)
        result = numpy.ndarray((self.nummet+1)*fT[0],dtype=numpy.float64)
        success = numpy.array(0,dtype=numpy.int32,ndmin=1)
        ret = self.sim_model(envar,
                             param,
                             fT,
                             control,
                             result,
                             success)
        ret = {
            'colT':result[0:fT[0]],
            'success':success
            }
        for n in range(self.nummet):
            ret[self.metnames[n]] = result[((n+1)*fT[0]):((n+2)*fT[0])]
        return ret

"""
TEST
"""
if __name__ == '__main__':
    pass