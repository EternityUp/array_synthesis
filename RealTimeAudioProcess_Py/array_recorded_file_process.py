# -*- encoding=utf-8 -*-

import wave
import numpy as np
import array_synthesis

print dir(array_synthesis)

# open wave file to be processed
wf = wave.open(u"./xmos_record1_chs4.wav", "rb")
wf_out = wave.open(u"./xmos_record1_ch_out.wav", "wb")

print dir(wf)

sample_rate = wf.getframerate()
num_channels = wf.getnchannels()
num_frames = wf.getnframes()
sample_width = wf.getsampwidth()
print "sample_rate={},num_channels={},num_frames={},sample_width={}"\
    .format(sample_rate, num_channels, num_frames, sample_width)

wf_out.setframerate(sample_rate)
wf_out.setnchannels(1)
wf_out.setsampwidth(sample_width)


ret = 0
frame_num = 0
single_ch_frame_len = 160
multi_chs_frame_len = single_ch_frame_len * num_channels


def Int16_ToFloat(in_int16_data, channels, channel_length):
    out_norm_float_data = np.zeros((channels, channel_length), dtype=np.float32)
    for i in range(channels):
        for j in range(channel_length):
            index = j * channels + i
            if in_int16_data[index] > 0:
                out_norm_float_data[i, j] = in_int16_data[index] / 32767.
            else:
                out_norm_float_data[i, j] = in_int16_data[index] / 32768.
    return out_norm_float_data


def Float32_To_Int16(float32_audio, length):
    norm_int16_audio = np.zeros(length, dtype=np.int16)
    for i in range(length):
        if float32_audio[i] > 0:
            norm_int16_audio[i] = float32_audio[i] * 32767
        else:
            norm_int16_audio[i] = float32_audio[i] * 32768
    return norm_int16_audio


py_aps_inst = array_synthesis.pyArrayProcessSynthesis()
print dir(py_aps_inst)
# 初始化
py_aps_inst.PyInitArrayProcessInst()
py_aps_inst.PyInitMicArray()
# 设置参数
py_aps_inst.PySetParas()
# 验证参数
py_aps_inst.PyValidateParas()

out_single_ch_audio = np.zeros(single_ch_frame_len, dtype=np.float32)
int16_out_single_ch_audio = np.zeros(single_ch_frame_len, dtype=np.int16)
while ret < num_frames:
    audio_str = wf.readframes(single_ch_frame_len)  # 读取字符串格式的音频数据
    read_samples = len(audio_str) / num_channels / sample_width
    int16_audio_data = np.fromstring(audio_str, dtype=np.int16)  # 将字符串转化成int16类型
    # 分离多通道音频数据并且进行数值归一化转化成float32类型
    in_float_audio_data = Int16_ToFloat(int16_audio_data, num_channels, read_samples)

    # 处理核心
    py_aps_inst.PyArrayProcessCore(in_float_audio_data, out_single_ch_audio)

    # 将数据从float32类型转化成int16类型
    int16_out_single_ch_audio = Float32_To_Int16(out_single_ch_audio, single_ch_frame_len)
    wf_out.writeframes(int16_out_single_ch_audio.tostring())

    ret += read_samples
    frame_num += 1


# 释放资源
py_aps_inst.PyFreeProcInst()



