clear all
close all
clc
% 根据xmos硬件多通道录音构造阵列信号处理输入文件
[audio1,fs] = audioread('./xmos录音/mic1.wav');
[audio2,~] = audioread('./xmos录音/mic3.wav');
[audio3,~] = audioread('./xmos录音/mic4.wav');
[audio4,~] = audioread('./xmos录音/mic6.wav');
audio = [audio1 audio2 audio3 audio4];
audiowrite('xmos_record1_chs4.wav',audio,fs);


[audio12,fs] = audioread('./xmos录音/mic1_2.wav');
[audio22,~] = audioread('./xmos录音/mic3_2.wav');
[audio32,~] = audioread('./xmos录音/mic4_2.wav');
[audio42,~] = audioread('./xmos录音/mic6_2.wav');
audio_2 = [audio12 audio22 audio32 audio42];
audiowrite('xmos_record2_chs4.wav',audio_2,fs);
