clear all
close all
clc
M = 2048;
x1 = randn(M,1);
x2 = randn(M,1); 
nfft = 256;
noverlap = 0;
%% ���������źŵ������
% ����׵ļ�����һ��ͳ�ƹ���
% ����Ϊ�������׺��Թ����׵ļ������ͳ�ƽ��
% �����֡�����֡����������ף����������
% ֡��Խ�࣬ͳ�ƽ��Խ�ɿ�������׽��Խ�ӽ�����ƫ����
mscohere(x1,x2,hanning(nfft),noverlap,nfft); % Plot estimate

