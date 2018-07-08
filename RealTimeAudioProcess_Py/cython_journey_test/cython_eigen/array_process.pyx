cimport numpy as np
import ctypes

cdef extern from "./inc/array_process_wrapper.h":
    cdef cppclass ArrayProcessSynthesis:
        ArrayProcessSynthesis()
        void InitArrayProcessSynthesis()
        void InitMicArray()
        void SetParas()
        void ValidateParas()
        void ArrayProcessCore(float *in_multi_chs_data, float *out_single_ch_data)
        void FreeProcInst()
        float **in_chs_audio_data
        float *out_ch_audio_data


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

  def PyArrayProcessCore(self, np.ndarray[np.float32_t, ndim=2] in_audio, np.ndarray[np.float32_t, ndim=1] out_audio):
      self.thisptr.ArrayProcessCore(<float*>in_audio.data, <float*>out_audio.data)

  def PyFreeProcInst(self):
      self.thisptr.FreeProcInst()

  def printHelloTest(self):
    print "hello, class pyArrayProcessSynthesis!"