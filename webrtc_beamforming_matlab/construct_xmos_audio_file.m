clear all
close all
clc
% ����xmosӲ����ͨ��¼�����������źŴ��������ļ�
[audio1,fs] = audioread('./xmos¼��/mic1.wav');
[audio2,~] = audioread('./xmos¼��/mic3.wav');
[audio3,~] = audioread('./xmos¼��/mic4.wav');
[audio4,~] = audioread('./xmos¼��/mic6.wav');
audio = [audio1 audio2 audio3 audio4];
audiowrite('xmos_record1_chs4.wav',audio,fs);


[audio12,fs] = audioread('./xmos¼��/mic1_2.wav');
[audio22,~] = audioread('./xmos¼��/mic3_2.wav');
[audio32,~] = audioread('./xmos¼��/mic4_2.wav');
[audio42,~] = audioread('./xmos¼��/mic6_2.wav');
audio_2 = [audio12 audio22 audio32 audio42];
audiowrite('xmos_record2_chs4.wav',audio_2,fs);
