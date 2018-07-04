#! /usr/bin/env python
# -*- encoding=utf-8 -*-

import wave
from ctypes import *
import os
import sys
import numpy as np

env_dist = os.environ


if 'LD_LIBRARY_PATH' not in os.environ:
    os.environ['LD_LIBRARY_PATH'] = os.getcwd() + u"/libs/"
    try:
        os.execv(sys.argv[0], sys.argv)
    except Exception, exc:
        print 'Failed re-exec:', exc
        sys.exit(1)

env_dist = os.environ
print env_dist.get('LD_LIBRARY_PATH')

lib = CDLL(os.getcwd() + u"/libs/libArrayProcessSynthesis.so")
lib2 = CDLL(os.getcwd() + u"/libs/libAudioFileProcess.so")

# open wave file to be processed
wf = wave.open(u"./xmos_record1_chs4.wav","rb")
print dir(wf)

sample_rate = wf.getframerate()
num_channels = wf.getnchannels()
num_frames = wf.getnframes()
sample_width = wf.getsampwidth()
print "sample_rate={},num_channels={},num_frames={},sample_width={}"\
    .format(sample_rate, num_channels, num_frames, sample_width)

ret = 0
frame_num = 0
single_ch_frame_len = 160
multi_chs_frame_len = single_ch_frame_len * num_channels

Int16_tToFloat = lib2.Int16_tToFloat
Int16_tToFloat.restype = None
process_core = lib.ProcCoreMicArrayProcSingleOut

# 定义各类待处理的结构体类型
class Mic_Array_Process_Inst(Structure):
    _fields_ = []











float_audio_data = np.asarray(np.zeros((multi_chs_frame_len, 1)), dtype=np.float32)
float_audio_data_ptr = cast(float_audio_data.ctypes.data, POINTER(c_float))

# 创建用于输出的1维float类型数组的指针
out_float_audio_data = np.asarray(np.zeros((single_ch_frame_len, 1)), dtype=np.float32)
out_float_audio_data_ptr = cast(out_float_audio_data.ctypes.data, POINTER(c_float))
# 创建用于输入的二维float类型的数组的指针
float_array_frame_len = (c_float * single_ch_frame_len)
float_2D_array_ptr = (POINTER(float_array_frame_len) * num_channels)
float_audio_2D_ptr = float_2D_array_ptr()  #

while ret <= num_frames:
    audio_str = wf.readframes(single_ch_frame_len)  # 读取字符串格式的音频数据
    read_samples = len(audio_str) / num_channels / sample_width
    int16_audio_data = np.fromstring(audio_str, dtype=np.int16)  # 将字符串转化成int16类型
    int16_audio_data_ptr = cast(int16_audio_data.ctypes.data, POINTER(c_short))  # 获取一维数组指针
    Int16_tToFloat(int16_audio_data_ptr, c_long(len(int16_audio_data)), float_audio_data_ptr)  # 数值归一化
    float_audio_data.shape = num_channels, -1  # 分离多通道数据, 转换成多维数组

    for i in range(num_channels):
        float_audio_2D_ptr[i] = pointer(float_array_frame_len(* float_audio_data[i].tolist()))

    ret += single_ch_frame_len
    frame_num += 1





