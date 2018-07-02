# -*- encoding = utf-8 -*-
from audio_recorder_wrapper import *
import time

record = Recorder(channels=6,rate=44100,frames_per_buffer=1024)

# non callback mode, namely blocking mode
with record.open("non_callback_record.wav", "wb") as record_file:
    record_file.record(duration=5.0)
    record_file.stop_recording_stream()


record_callback = Recorder(channels=6,rate=44100,frames_per_buffer=1024)

# callback mode, namely non blocking mode
with record_callback.open("callback_record.wav","wb") as record_file_callback:
    record_file_callback.record_in_callback_mode()
    record_file_callback.start_recording_stream()
    time.sleep(5.0)
    record_file_callback.stop_recording_stream()

