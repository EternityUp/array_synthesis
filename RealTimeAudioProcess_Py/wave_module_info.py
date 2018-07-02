# -*- encoding=utf-8 -*-
import wave
print dir(wave)

wf = wave.open("327_16.wav", mode="rb")
print dir(wf)

frame_rate = wf.getframerate()
channels = wf.getnchannels()
frames = wf.getnframes()
sample_width = wf.getsampwidth()

print "frame_rate={}, channels={}, frames={}, sample_width={}".format(frame_rate, channels, frames, sample_width)


