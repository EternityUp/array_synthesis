clear all
close all
clc
M = 2048;
x1 = randn(M,1);
x2 = randn(M,1); 
nfft = 256;
noverlap = 0;
%% 计算两个信号的相干谱
% 相干谱的计算是一个统计过程
% 体现为互功率谱和自功率谱的计算均是统计结果
% 如果单帧或较少帧数计算相干谱，会产生错误
% 帧数越多，统计结果越可靠，相干谱结果越接近于无偏估计
mscohere(x1,x2,hanning(nfft),noverlap,nfft); % Plot estimate

