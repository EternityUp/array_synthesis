from test import pyTest
import ctypes
import numpy as np
square = ctypes.cdll.LoadLibrary("./square.so")
n = 5
a = np.arange(n, dtype=np.double)
aptr = a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
square.square(aptr, n)
print a

k = pyTest(10)
k.returnFive()
print k + k
print k - k
print k.printMe()

m = k.getNDArray(4,4)
print m

mm = k.getNDArrayFromPtr(4,4)
print mm





