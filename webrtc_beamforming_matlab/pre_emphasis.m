clear all
close all
clc
[audio,fs] = audioread('stereo.wav');
audio1 = audio(:,1);
x = audio1(1:160);
alpha = 0.9375;
win = [1 -alpha];
y1 = filter(win,1,x);
y2 = [x(1);x(2:end) - alpha * x(1:end-1)];
diffy = y2 - y1;
max(abs(diffy))


x1 = [zeros(224,1);audio(1:32)'];
x2 = filter(win,1,x1);
x22 = [x1(1);x1(2:end) - alpha * x1(1:end-1)];