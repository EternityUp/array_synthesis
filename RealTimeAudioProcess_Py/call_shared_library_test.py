#! /usr/bin/env python
# -*- encoding=utf-8 -*-

from ctypes import *
import os
import sys
import numpy as np

env_dist = os.environ
print type(env_dist)
print env_dist

# 打印所有环境变量，遍历字典
for key in env_dist:
    print key + ' : ' + env_dist[key]

if 'LD_LIBRARY_PATH' not in os.environ:
    os.environ['LD_LIBRARY_PATH'] = os.getcwd() + u"/libs/"
    try:
        os.execv(sys.argv[0], sys.argv)
    except Exception, exc:
        print 'Failed re-exec:', exc
        sys.exit(1)

print env_dist.get('LD_LIBRARY_PATH')
print env_dist['LD_LIBRARY_PATH']


print os.getcwd()

lib1 = CDLL(u"./libs/libArrayShapeAnalysis.so")
print lib1, type(lib1)

lib2 = CDLL(os.getcwd() + u"/libs/libAudioFileProcess.so")
print lib2, type(lib2)

lib3 = CDLL(os.getcwd() + u"/libs/libBasicSignalProcess.so")

lib4 = CDLL(os.getcwd() + u"/libs/libFourierTransform.so")

lib5 = CDLL(os.getcwd() + u"/libs/libWavFileRW.so")

lib6 = CDLL(os.getcwd() + u"/libs/libArrayProcessSynthesis.so")


class Point(Structure):
    _fields_ = [("x", c_float),
                ("y", c_float),
                ("z", c_float)]


a = Point()
print a.x, a.y, a.z

a = Point(1.0, 2.0, 3.0)
print a.x, a.y, a.z

b = Point(1.0, 2.0, 3.0)

pa = pointer(a)
pb = pointer(b)

print pa, pb, pa.contents, pb.contents

DotProduct = lib1.DotProduct
DotProduct.restype = c_float
c = DotProduct(byref(a), byref(b))
print c

ComputeSoundPath = lib1.ComputeSoundPath
ComputeSoundPath.restype = c_float
d = ComputeSoundPath(byref(a), c_float(0.5))
print d


Int16_tToFloat = lib2.Int16_tToFloat
Int16_tToFloat.restype = None

src = (c_short * 2)(1,2)
print src[0], src[1]

dst = (c_float * 2)()
print dst[0], dst[1]

Int16_tToFloat(byref(src), c_long(len(src)), byref(dst))

print c_long(len(src)).value
print dst[0], dst[1]

print lib5, type(lib5)
print lib6, type(lib6)

WriteWavHeader = lib5.WriteWavHeader
print type(WriteWavHeader), WriteWavHeader

NormComplex = lib6.ComputeOutSpatialSpectrumMvdr
print type(NormComplex), NormComplex








