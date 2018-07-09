clear all
close all
clc
[audio1,fs] = audioread('stereo2_out1次加窗.wav');
[audio2,~] = audioread('stereo2_out2次加窗.wav');
[M,N] = size(audio1);
t = (1:M)/fs;

figure
hold on
plot(t,audio1(:,1),'b',t,audio2(:,1),'r');
