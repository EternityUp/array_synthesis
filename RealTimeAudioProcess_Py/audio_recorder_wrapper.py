# -*- coding: utf-8 -*-
import pyaudio
import wave


class Recorder(object):
    """
    a recorder class for recording audio to a wav file
    """

    def __init__(self, channels=2, rate=44100, frames_per_buffer=1024):
        self.channels = channels
        self.rate = rate
        self.frames_per_buffer = frames_per_buffer

    def open(self, filename, mode="wb"):
        return RecordingFile(filename, mode, self.channels, self.rate, self.frames_per_buffer)


class RecordingFile(object):
    def __init__(self, filename, mode, channels,
                 rate, frames_per_buffer):
        self.filename = filename
        self.mode = mode
        self.channels = channels
        self.rate = rate
        self.frames_per_buffer = frames_per_buffer

        self._pa = pyaudio.PyAudio()
        self.wav_file = self._prepare_file(self.filename, self.mode)
        self._stream = None

    def __enter__(self):
        return self

    def __exit__(self, exception, value, traceback):
        self.close()

    def record(self, duration):
        print "open audio stream in blocking mode"
        self._stream = self._pa.open(format=pyaudio.paInt16,
                                     channels=self.channels,
                                     rate=self.rate,
                                     input=True,
                                     frames_per_buffer=self.frames_per_buffer)
        print "start recording in blocking record ... "

        frame_num = int(self.rate / self.frames_per_buffer * duration)
        for i in range(frame_num):
            audio = self._stream.read(self.frames_per_buffer)
            self.wav_file.writeframes(audio)
        return None

    def record_in_callback_mode(self):
        print "open audio stream in non blocking mode "
        self._stream = self._pa.open(format=pyaudio.paInt16,
                                     channels=self.channels,
                                     rate=self.rate,
                                     input=True,
                                     frames_per_buffer=self.frames_per_buffer,
                                     stream_callback=self.write_to_wav_callback)

    def start_recording_stream(self):
        print "start audio stream"
        self._stream.start_stream()
        return self

    def stop_recording_stream(self):
        print "stop audio stream"
        self._stream.stop_stream()
        return self

    def write_to_wav_callback(self, in_data, frame_count, time_info, status):
        self.wav_file.writeframes(in_data)
        return in_data, pyaudio.paContinue

    def _prepare_file(self, filename, mode):
        wav_file = wave.open(filename, mode)
        wav_file.setnchannels(self.channels)
        wav_file.setframerate(self.rate)
        wav_file.setsampwidth(self._pa.get_sample_size(pyaudio.paInt16))
        return wav_file

    def close(self):
        self._stream.close()
        self._pa.terminate()
        self.wav_file.close()




