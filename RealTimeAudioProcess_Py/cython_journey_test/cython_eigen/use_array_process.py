import array_synthesis as a_s
import numpy as np
import sys

print dir(a_s)

x = a_s.pyArrayProcessSynthesis()

print x

x.printHelloTest()


print dir(x)


x.PyInitArrayProcessInst()

x.PySetParas()

x.PyValidateParas()

x.PyInitMicArray()


in_multi_chs_audio = np.random.random((4, 160))
print in_multi_chs_audio.dtype

out_single_ch_audio = np.zeros((1, 160), dtype=np.float32)
print out_single_ch_audio.dtype


x.PyArrayProcessCore(in_multi_chs_audio, out_single_ch_audio)


x.PyFreeProcInst()




