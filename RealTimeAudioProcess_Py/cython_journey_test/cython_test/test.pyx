include "matrixxfpy.pyx"
import numpy

cdef extern from "cpp_test.h":
  cdef cppclass Test:
    Test()
    Test(int test1)
    int test1
    int returnFive()
    Test add "operator+"(Test other)
    Test sub "operator-"(Test other)
    MatrixXfPy *ptr;
    MatrixXfPy getMatrixXf(int d1, int d2)
    MatrixXfPy* getMatrixXfFromPtr(MatrixXfPy *ptr, int d1, int d2);
    
cdef class pyTest:
  cdef Test* thisptr
  def __cinit__(self, int test1):
    self.thisptr = new Test(test1)
  def __dealloc__(self):
    del self.thisptr
  
  def __add__(pyTest left, pyTest other):
    cdef Test t = left.thisptr.add(other.thisptr[0])
    cdef pyTest tt = pyTest(t.test1)
    return tt
  
  
  def __sub__(pyTest left, pyTest other):
    cdef Test t = left.thisptr.sub(other.thisptr[0])
    cdef pyTest tt = pyTest(t.test1)
    return tt
  
  
  def __repr__(self):
    return "pyTest[%s]" % (self.thisptr.test1)
  
  def returnFive(self):
    return self.thisptr.returnFive()
  
  def printMe(self):
    return "hello world"
  
  def getNDArray(self, int d1, int d2): 
    cdef MatrixXfPy me = self.thisptr.getMatrixXf(d1,d2)  
    result = numpy.zeros((me.rows(),me.cols())) 
    for row in range(me.rows()): 
        for col in range(me.cols()): 
            result[row, col] = me.coeff(row, col)   
    return result 
  
  def getNDArrayFromPtr(self, int d1, int d2):
    self.thisptr.ptr = self.thisptr.getMatrixXfFromPtr(self.thisptr.ptr, d1, d2)
    print "hello, now in getNDArrayFromPtr function!"
    print self.thisptr.ptr.rows(), self.thisptr.ptr.cols()
    
    result = numpy.zeros((self.thisptr.ptr.rows(),self.thisptr.ptr.cols())) 
    for row in range(self.thisptr.ptr.rows()): 
        for col in range(self.thisptr.ptr.cols()): 
            result[row, col] = self.thisptr.ptr.coeff(row, col)   
    return result 
  
  
  
  
  
  
  
  
  
  