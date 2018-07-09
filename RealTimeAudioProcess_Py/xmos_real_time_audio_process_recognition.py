# -*- encoding=utf-8 -*-

import sys
import wave
import pyaudio
import speech_recognition as sr
import threading
import array_synthesis
import numpy as np

sample_frequency = 16000
channels = 6  # 2,3,4,5通道为原始音频流,通道1为xmos输出音频流，通道6为回声参考通道
frame_len = 160
duration = 10
sample_width = 2

# 全局变量定义
raw_audio_data = None
processed_audio_data = None
audio_record_prepared = False
audio_process_prepared = False
original_frame = []
processed_frame = []
xmos_out_frame = []

raw_audio_file = wave.open(sys.argv[1], "wb")
raw_audio_file.setsampwidth(sample_width)
raw_audio_file.setframerate(sample_frequency)
raw_audio_file.setnchannels(channels)

# 录制一定时长的实时音频，并且写入.wav文件
def audio_record(sample_rate, num_channels, chunk_size, seconds):
    global raw_audio_data, audio_record_prepared, original_frame
    print "Thread (name:%s) for recording audio is running ..." % threading.current_thread().name
    pa = pyaudio.PyAudio()
    print "open audio stream for recording ... "
    stream = pa.open(format=pyaudio.paInt16,
                     channels=num_channels,
                     rate=sample_rate,
                     input=True,
                     frames_per_buffer=chunk_size)
    print "open a wav file for writing original audio stream"
    frame_num = sample_rate / chunk_size * seconds
    print "recording ... "
    for i in range(frame_num):
        raw_audio_data = stream.read(chunk_size)
        original_frame.append(raw_audio_data)
        raw_audio_file.writeframes(raw_audio_data)
    audio_record_prepared = True  # 待处理音频准备完毕，标志位置为True
    stream.stop_stream()
    stream.close()
    pa.terminate()


audio_process_inst = array_synthesis.pyArrayProcessSynthesis()


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


def Init_Audio_Process():
    # 初始化
    audio_process_inst.PyInitArrayProcessInst()
    audio_process_inst.PyInitMicArray()
    # 设置参数
    audio_process_inst.PySetParas()
    # 验证参数
    audio_process_inst.PyValidateParas()


processed_audio_file = wave.open(sys.argv[2], "wb")
processed_audio_file.setframerate(sample_frequency)
processed_audio_file.setnchannels(1)
processed_audio_file.setsampwidth(sample_width)

xmos_out_file = wave.open(sys.argv[3], "wb")
xmos_out_file.setframerate(sample_frequency)
xmos_out_file.setnchannels(1)
xmos_out_file.setsampwidth(sample_width)



# 处理已经录制完毕的音频流
def process_audio(num_channels, frame_length):
    global raw_audio_data, processed_audio_data, audio_record_prepared, audio_process_prepared, processed_frame
    print "Thread (name:%s) for processing audio is running ..." % threading.current_thread().name
    while True:
        if audio_record_prepared:
            print "start processing original audio stream ... "
            in_audio_stream = original_frame
            frame_num = len(in_audio_stream)  # 待处理的音频流帧数
            # 分帧处理
            out_single_ch_audio = np.zeros(frame_length, dtype=np.float32)
            Init_Audio_Process()
            for i in range(frame_num):
                frame_audio_bytes = in_audio_stream[i]
                int16_audio_data = np.fromstring(frame_audio_bytes, dtype=np.int16)  # 将字符串转化成int16类型
                # 分离多通道音频数据并且进行数值归一化转化成float32类型
                in_float_audio_data = Int16_ToFloat(int16_audio_data, num_channels, frame_length)
                in_float_audio_data_m4 = in_float_audio_data[1:5, :]
                # 处理核心
                audio_process_inst.PyArrayProcessCore(in_float_audio_data_m4, out_single_ch_audio)
                # 将数据从float32类型转化成int16类型
                int16_out_single_ch_audio = Float32_To_Int16(out_single_ch_audio, frame_length)
                # 施加一定的幅度增益
                int16_out_single_ch_audio = 10 * int16_out_single_ch_audio
                processed_audio_file.writeframes(int16_out_single_ch_audio.tostring())
                processed_frame.append(int16_out_single_ch_audio.tostring())
                # 提取第1通道的音频
                xmos_out_ch_data = Float32_To_Int16(in_float_audio_data[0, :], frame_length)
                xmos_out_frame.append(xmos_out_ch_data.tostring())
                xmos_out_file.writeframes(xmos_out_ch_data.tostring())
            audio_process_prepared = True
            break


# 音频流处理完毕后，基于谷歌和微软识别引擎识别语音
def recognize_processed_audio(sample_rate, bytes_per_frame):
    global processed_frame
    print "Thread (name:%s) for recognizing audio is running ..." % threading.current_thread().name
    while True:
        if audio_process_prepared:
            r = sr.Recognizer()
            audio = sr.AudioData(b''.join(processed_frame), sample_rate, bytes_per_frame)
            print "start recognizing processed audio ... "
            # recognize speech using Google Speech Recognition  测试谷歌语音识别引擎
            try:
                # for testing purpose,we're just using the default API key
                # to use another API key, use 'r.recognize_google(audio,key="GOOGLE_SPEECH_RECOGNITION_API_KEY")'
                # instead of 'r.recognize_google(audio)
                print("Google Speech Recognition thinks you said " + r.recognize_google(audio))
            except sr.UnknownValueError:
                print("Google Speech Recognition could not understand audio")
            except sr.RequestError as e:
                print("Could not request results from Google Speech Recognition service;{0}".format(e))

                # recognize speech using Microsoft Bing Voice Recognition
                # Microsoft Bing Voice Recognition API keys 32-character lowercase hexadecimal strings
            BING_KEY = "6c0b8b2e8bcd4b138c97ad103d6001ca"
            try:
                print("Microsoft Bing Voice Recognition thinks you said " + r.recognize_bing(audio, key=BING_KEY))
            except sr.UnknownValueError:
                print("Microsoft Bing Voice Recognition could not understand audio")
            except sr.RequestError as e:
                print("Could not request results from Microsoft Bing Voice Recognition service; {0}".format(e))

            # recognize speech using Google Cloud Speech
            GOOGLE_CLOUD_SPEECH_CREDENTIALS = r"""{
              "type": "service_account",
              "project_id": "able-balm-186710",
              "private_key_id": "23a05ee9d71d1539e93bc47e11190898488263a8",
              "private_key": "-----BEGIN PRIVATE KEY-----\nMIIEvAIBADANBgkqhkiG9w0BAQEFAASCBKYwggSiAgEAAoIBAQCEvB7RxAdA8F3Z\neCTbm8BEpux0eLMa4APyxukxkHbJZ9bxVHswkkkLVpWaQmXYTJU1pK00UilJ2Eou\nyEgG/q0mcV1KflQw0OMNZKty7nDrRRVhWxqrS9x2azA7rN3m/WFJ+ZnLvcRZRSFS\nMEDwMCtmRJnp7UTq1HbzFsE7AJEqbAaF260VWDSpoKNHZaERIK5Jd4pQRN7WhbOz\ne2ItznvENszxC8eY4k5Wx3uuB21EjV9yYrmp0K/yw9npA+Wp14+8NymbCRwLWb8o\nMI09dudvS9GSpzl0zrrGHUFx9I8usHUrEIXu/pSUf3+uy8t8O9SbE0xumhq9W7SB\nt4iN88JpAgMBAAECggEANdd5SO2TY6wc5ZHuoYZUP4TWdYhgXo5wC5LUFN7c0GXm\nso3qSpGMP6UnmiNE3LBz9gaOm6KYTPQHY2KqlOvJoMZraURFMmgJMe2b/krqUwT6\n3zrtXn6vNvnT3uEIQiKVzEQLNrILa4q3MXeRX9yqPhmltVYhloAIoCKizyQzvljt\n39Fk6vF8RAvJdrEYbI5AlK8BFppvzxkuUfUPdaK821OB+zrhRKbNjqCo/09d6u+/\nVHLG6TUDFP4KM9DhgC6HHWn/sS8d8PVGhOJB6U/tCij1H6X96PnQDIYHy5CIYxvu\npGe+gvoDf1LEy4QDaJYyccrjECl856tS2SMSKNqbYQKBgQC40HN+Fg/4dbCEaISr\niyuwiOUMhLaDXhW3oy+aoyaVvNoUvHtxQ8zRzbHYZzcHNu/Vnz0zL1WcJlc72PN3\n3cnrCjpF7aEtBjyZYbLyeBeSe5VrTI/MoK7ovjx5t2YWfVZQVmDhPrmjQevccQSI\n+jlay3/4dVTMOSVfGaVTTw7e2wKBgQC33GVJnGE0D3+DiVdMmYmgbhFRmQ/MPZH7\nzWTNrcOvgQDC2U2UC4LcPhDdcPNmpbZmbG787P/tkWbY5wzBwlib3/gkDRmEPlKe\nFtWr6pB/nIYaYn0jzkWi+eAfTnfIs3e1vNHSCciKJ5sSHeBigv7nIQf+vUprZRYb\nbi1e9oI9CwKBgCA8Q5i/cmOs62/85v8g1CsDhagMUeVR+MnCDeCUCnvdp9AsG//i\niowhq56KSj/Y5jZFgyA1ZmJJEbgfnD/REJINg6KE0zMOPm2ma9b8+WeUZLiFbyOK\n91cjL1svkP/lNrPmjRlcnnaKXgBiOh9GPdDHY/fLR1IjuY//4iVIydg3AoGANBMm\nVP0HwvkIwraplTZ/doL2QMg7YzulF06LWH34yMKe9pEZme7Qt5SUrkOJjO8uhD/+\nB4EQ07a1DIYUZOSouC1tWLilG4GTYNdS2YhsONiaWuq+St/ndUuUoQlWf+/k3gKG\n4xiGRisFjtILdZtomwoN+6adZ2GQK2C/VQA4zxUCgYA+XF2xcWXrbudjpzSNO1tn\nT8mJWWWj8W15t/nDSukifo1KNo8LOTf2s1doDZJbo/78mNZoZFAiGnLcrnUtVKTO\n+d0p2FatQs4rvzPRb5s76LssyVRfK6mTMdIGc6f6+23OPh/PyEo3BZ6c9a9H03vu\ngOOIRw8SR457WzL2ByoyVw==\n-----END PRIVATE KEY-----\n",
              "client_email": "amy-speech@able-balm-186710.iam.gserviceaccount.com",
              "client_id": "117483853940217018962",
              "auth_uri": "https://accounts.google.com/o/oauth2/auth",
              "token_uri": "https://accounts.google.com/o/oauth2/token",
              "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
              "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/amy-speech%40able-balm-186710.iam.gserviceaccount.com"
            }
            """
            try:
                print("Google Cloud Speech thinks you said " + r.recognize_google_cloud(audio,
                                                                    credentials_json=GOOGLE_CLOUD_SPEECH_CREDENTIALS))
            except sr.UnknownValueError:
                print("Google Cloud Speech could not understand audio")
            except sr.RequestError as e:
                print("Could not request results from Google Cloud Speech Service;{0}".format(e))
            break


# 基于谷歌和微软识别引擎识别xmos阵列算法处理输出的语音
def recognize_xmos_out_audio(sample_rate, bytes_per_frame):
    global xmos_out_frame
    print "Thread (name:%s) for recognizing audio is running ..." % threading.current_thread().name
    while True:
        if audio_process_prepared:
            r = sr.Recognizer()
            audio = sr.AudioData(b''.join(xmos_out_frame), sample_rate, bytes_per_frame)
            print "(for XMOS)start recognizing processed audio ... "
            # recognize speech using Google Speech Recognition  测试谷歌语音识别引擎
            try:
                # for testing purpose,we're just using the default API key
                # to use another API key, use 'r.recognize_google(audio,key="GOOGLE_SPEECH_RECOGNITION_API_KEY")'
                # instead of 'r.recognize_google(audio)
                print("(for XMOS)Google Speech Recognition thinks you said " + r.recognize_google(audio))
            except sr.UnknownValueError:
                print("(for XMOS)Google Speech Recognition could not understand audio")
            except sr.RequestError as e:
                print("(for XMOS)Could not request results from Google Speech Recognition service;{0}".format(e))

                # recognize speech using Microsoft Bing Voice Recognition
                # Microsoft Bing Voice Recognition API keys 32-character lowercase hexadecimal strings
            BING_KEY = "6c0b8b2e8bcd4b138c97ad103d6001ca"
            try:
                print("(for XMOS)Microsoft Bing Voice Recognition thinks you said " + r.recognize_bing(audio, key=BING_KEY))
            except sr.UnknownValueError:
                print("(for XMOS)Microsoft Bing Voice Recognition could not understand audio")
            except sr.RequestError as e:
                print("(for XMOS)Could not request results from Microsoft Bing Voice Recognition service; {0}".format(e))

            # recognize speech using Google Cloud Speech
            GOOGLE_CLOUD_SPEECH_CREDENTIALS = r"""{
              "type": "service_account",
              "project_id": "able-balm-186710",
              "private_key_id": "23a05ee9d71d1539e93bc47e11190898488263a8",
              "private_key": "-----BEGIN PRIVATE KEY-----\nMIIEvAIBADANBgkqhkiG9w0BAQEFAASCBKYwggSiAgEAAoIBAQCEvB7RxAdA8F3Z\neCTbm8BEpux0eLMa4APyxukxkHbJZ9bxVHswkkkLVpWaQmXYTJU1pK00UilJ2Eou\nyEgG/q0mcV1KflQw0OMNZKty7nDrRRVhWxqrS9x2azA7rN3m/WFJ+ZnLvcRZRSFS\nMEDwMCtmRJnp7UTq1HbzFsE7AJEqbAaF260VWDSpoKNHZaERIK5Jd4pQRN7WhbOz\ne2ItznvENszxC8eY4k5Wx3uuB21EjV9yYrmp0K/yw9npA+Wp14+8NymbCRwLWb8o\nMI09dudvS9GSpzl0zrrGHUFx9I8usHUrEIXu/pSUf3+uy8t8O9SbE0xumhq9W7SB\nt4iN88JpAgMBAAECggEANdd5SO2TY6wc5ZHuoYZUP4TWdYhgXo5wC5LUFN7c0GXm\nso3qSpGMP6UnmiNE3LBz9gaOm6KYTPQHY2KqlOvJoMZraURFMmgJMe2b/krqUwT6\n3zrtXn6vNvnT3uEIQiKVzEQLNrILa4q3MXeRX9yqPhmltVYhloAIoCKizyQzvljt\n39Fk6vF8RAvJdrEYbI5AlK8BFppvzxkuUfUPdaK821OB+zrhRKbNjqCo/09d6u+/\nVHLG6TUDFP4KM9DhgC6HHWn/sS8d8PVGhOJB6U/tCij1H6X96PnQDIYHy5CIYxvu\npGe+gvoDf1LEy4QDaJYyccrjECl856tS2SMSKNqbYQKBgQC40HN+Fg/4dbCEaISr\niyuwiOUMhLaDXhW3oy+aoyaVvNoUvHtxQ8zRzbHYZzcHNu/Vnz0zL1WcJlc72PN3\n3cnrCjpF7aEtBjyZYbLyeBeSe5VrTI/MoK7ovjx5t2YWfVZQVmDhPrmjQevccQSI\n+jlay3/4dVTMOSVfGaVTTw7e2wKBgQC33GVJnGE0D3+DiVdMmYmgbhFRmQ/MPZH7\nzWTNrcOvgQDC2U2UC4LcPhDdcPNmpbZmbG787P/tkWbY5wzBwlib3/gkDRmEPlKe\nFtWr6pB/nIYaYn0jzkWi+eAfTnfIs3e1vNHSCciKJ5sSHeBigv7nIQf+vUprZRYb\nbi1e9oI9CwKBgCA8Q5i/cmOs62/85v8g1CsDhagMUeVR+MnCDeCUCnvdp9AsG//i\niowhq56KSj/Y5jZFgyA1ZmJJEbgfnD/REJINg6KE0zMOPm2ma9b8+WeUZLiFbyOK\n91cjL1svkP/lNrPmjRlcnnaKXgBiOh9GPdDHY/fLR1IjuY//4iVIydg3AoGANBMm\nVP0HwvkIwraplTZ/doL2QMg7YzulF06LWH34yMKe9pEZme7Qt5SUrkOJjO8uhD/+\nB4EQ07a1DIYUZOSouC1tWLilG4GTYNdS2YhsONiaWuq+St/ndUuUoQlWf+/k3gKG\n4xiGRisFjtILdZtomwoN+6adZ2GQK2C/VQA4zxUCgYA+XF2xcWXrbudjpzSNO1tn\nT8mJWWWj8W15t/nDSukifo1KNo8LOTf2s1doDZJbo/78mNZoZFAiGnLcrnUtVKTO\n+d0p2FatQs4rvzPRb5s76LssyVRfK6mTMdIGc6f6+23OPh/PyEo3BZ6c9a9H03vu\ngOOIRw8SR457WzL2ByoyVw==\n-----END PRIVATE KEY-----\n",
              "client_email": "amy-speech@able-balm-186710.iam.gserviceaccount.com",
              "client_id": "117483853940217018962",
              "auth_uri": "https://accounts.google.com/o/oauth2/auth",
              "token_uri": "https://accounts.google.com/o/oauth2/token",
              "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
              "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/amy-speech%40able-balm-186710.iam.gserviceaccount.com"
            }
            """
            try:
                print("(for XMOS)Google Cloud Speech thinks you said " + r.recognize_google_cloud(audio,
                                                                    credentials_json=GOOGLE_CLOUD_SPEECH_CREDENTIALS))
            except sr.UnknownValueError:
                print("(for XMOS)Google Cloud Speech could not understand audio")
            except sr.RequestError as e:
                print("(for XMOS)Could not request results from Google Cloud Speech Service;{0}".format(e))
            break


# create a thread for audio recording
audio_record_params = (sample_frequency, channels, frame_len, duration)
print audio_record_params
audio_record_thread = threading.Thread(target=audio_record, name="audio_record_thread", args=audio_record_params)

# create a thread for processing audio via mic array processing algorithms
audio_process_params = (channels, frame_len)
process_audio_thread = threading.Thread(target=process_audio, name="audio_process_thread", args=audio_process_params)


# create a thread for processed recorded audio recognition
recognize_audio_params = (sample_frequency, sample_width)
recognize_processed_audio_thread = threading.Thread(target=recognize_processed_audio, name="audio_recognize_thread",
                                          args=recognize_audio_params)

# create a thread for processed recorded audio recognition
recognize_xmos_out_audio_thread = threading.Thread(target=recognize_xmos_out_audio, name="audio_recognize_xmos_thread",
                                          args=recognize_audio_params)


# start and join created thread
audio_record_thread.start()
process_audio_thread.start()
recognize_processed_audio_thread.start()
recognize_xmos_out_audio_thread.start()
audio_record_thread.join()
process_audio_thread.join()
recognize_processed_audio_thread.join()
recognize_xmos_out_audio_thread.join()


raw_audio_file.close()
processed_audio_file.close()
xmos_out_file.close()











