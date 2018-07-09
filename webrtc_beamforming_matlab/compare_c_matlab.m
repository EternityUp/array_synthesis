% 比较c++和matlab两种语言实现的算法的差异
clear all
close all
clc
[audio_matlab,fs] = audioread('stereo2_out.wav');
[audio_c,~] = audioread('out.wav');
L = length(audio_c);
t = ( 1 : L ) / fs;
figure
subplot 211
hold on
plot(t,audio_matlab(:,1),'b.')
plot(t,audio_c(:,1),'r')
xlim([min(t) max(t)])
legend('matlab implementation','c++ implementation')
title('通道1')
subplot 212
hold on
plot(t,audio_matlab(:,2),'b.')
plot(t,audio_c(:,2),'r')
xlim([min(t) max(t)])
legend('matlab implementation','c++ implementation')
title('通道2')



figure
subplot 211
hold on
plot(t,audio_matlab(:,1),'b.')
plot(t,audio_c(:,1),'r')
xlim([4.1 4.15])
legend('matlab implementation','c++ implementation')
title('通道1局部放大对比图')
subplot 212
hold on
plot(t,audio_matlab(:,2),'b.')
plot(t,audio_c(:,2),'r')
xlim([4.1 4.15])
legend('matlab implementation','c++ implementation')
title('通道2局部放大对比图')