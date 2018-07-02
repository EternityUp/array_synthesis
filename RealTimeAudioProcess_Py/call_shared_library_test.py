# -*- encoding=utf-8 -*-
import ctypes
import os

print os.getcwd()

lib1 = ctypes.CDLL(u"./libs/libArrayShapeAnalysis.so")
lib2 = ctypes.CDLL(u"./libs/libAudioFileProcess.so")
lib3 = ctypes.CDLL(u"./libs/libBasicSignalProcess.so")
lib4 = ctypes.CDLL(u"./libs/libFourierTransform.so")
lib5 = ctypes.CDLL(u"./libs/libWavFileRW.so")

print type(lib1)

print dir(lib1)

print lib1.__dict__, lib1.__doc__
print lib1.__module__





