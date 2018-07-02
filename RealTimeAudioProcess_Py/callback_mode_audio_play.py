# -*- encoding=utf-8 -*-
import pyaudio
import wave
import time
import sys

if len(sys.argv) < 2:
    print("Play a wave file.\n\nUsage: %s filename.wav" % sys.argv[0])
    sys.exit(-1)

wf = wave.open(sys.argv[1],"rb")

p = pyaudio.PyAudio()


# define callback
def callback(in_data, frame_count, time_info, status):
    data = wf.readframes(frame_count)
    return data, pyaudio.paContinue


# open stream using callback (3)
stream = p.open(format=p.get_format_from_width(wf.getsampwidth()),
                channels=wf.getnchannels(),
                rate=wf.getframerate(),
                output=True,
                stream_callback=callback)

print dir(stream)

# start the stream (4)
stream.start_stream()

# wait for stream to finish (5)
while stream.is_active():
    print "stream.is_active(): %s" % stream.is_active()
    print "stream.is_stopped(): %s" % stream.is_stopped()
    time.sleep(0.1)


print "stream.is_active(): %s" % stream.is_active()
print "stream.is_stopped(): %s" % stream.is_stopped()

# stop stream (6)
stream.stop_stream()

print "stream.is_stopped(): %s" % stream.is_stopped()

stream.close()
wf.close()


# close PyAudio (7)
p.terminate()