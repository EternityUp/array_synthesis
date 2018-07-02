# -*- encoding=utf-8 -*-

import pyaudio

print dir(pyaudio)

print pyaudio.get_portaudio_version()
print pyaudio.get_portaudio_version_text()


audio_hdl = pyaudio.PyAudio()

print dir(audio_hdl)


default_host_api_info = audio_hdl.get_default_host_api_info()
print default_host_api_info

default_input_device_info = audio_hdl.get_default_input_device_info()
print default_input_device_info

default_output_device_info = audio_hdl.get_default_output_device_info()
print default_output_device_info

device_count = audio_hdl.get_device_count()
print device_count

host_api_count = audio_hdl.get_host_api_count()
print host_api_count

host_api_info_by_index = audio_hdl.get_host_api_info_by_index(1)
print host_api_info_by_index

device_info_by_index = audio_hdl.get_device_info_by_index(0)
print device_info_by_index

device_info_by_host_api_device_index = audio_hdl.get_device_info_by_host_api_device_index(0,0)
print device_info_by_host_api_device_index

host_api_info_by_type = audio_hdl.get_host_api_info_by_type(8)
print host_api_info_by_type

sample_size = audio_hdl.get_sample_size(8)
print "sample_size=%d" % sample_size

format_from_width = audio_hdl.get_format_from_width(2,unsigned=False)
print format_from_width

format_supported = audio_hdl.is_format_supported(rate=44100,input_device=0,input_channels=2,input_format=8,
                                                 output_device=0,output_channels=2,output_format=8)
print format_supported


