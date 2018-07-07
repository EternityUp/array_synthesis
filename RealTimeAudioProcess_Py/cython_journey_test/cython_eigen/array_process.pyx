import numpy as np

cdef extern from "./inc/array_process_wrapper.h":
    cdef cppclass ArrayProcessSynthesis:
        ArrayProcessSynthesis()
        void InitArrayProcessSynthesis()
        void InitMicArray()
        void SetParas()
        void ValidateParas()
        void ArrayProcessCore(float **in_multi_chs_data, float *out_single_ch_data)
        void FreeProcInst();


cdef class pyArrayProcessSynthesis:
  cdef ArrayProcessSynthesis* thisptr
  def __cinit__(self):
    self.thisptr = new ArrayProcessSynthesis()
  def __dealloc__(self):
    del self.thisptr

  def PyInitArrayProcessInst(self):
      self.thisptr.InitArrayProcessSynthesis()

  def PyInitMicArray(self):
      self.thisptr.InitMicArray()

  def PySetParas(self):
      self.thisptr.SetParas()

  def PyValidateParas(self):
      self.thisptr.ValidateParas()

  def PyArrayProcessCore(self, in_multi_chs_audio, out_single_ch_audio):
      cdef float **in_multi_chs_data
      cdef float *out_single_ch_data
      self.thisptr.ArrayProcessCore(in_multi_chs_data, out_single_ch_data)

  def PyFreeProcInst(self):
      self.thisptr.FreeProcInst()

  def printHelloTest(self):
    print "hello, class pyArrayProcessSynthesis!"