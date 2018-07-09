clear all
close all
clc
N = 1024;
x1 = randn( N ,1 );
M = 256;
delay = 20;
x11 = x1( 1 : M );
x12 = x1( 1 + delay : M + delay );
%% ʱ����ƻ����  �����ӳ���
[xc1,~] = xcorr(x11, x12);
[xc2,lags] = xcorr(x12, x11);
figure
subplot 211
plot(lags,xc1)
axis tight
subplot 212
plot(lags,xc2)
axis tight

[~,maxind1] = max(xc1);
delay_est1 = maxind1 - M; 

[maxval,maxind2] = max(xc2);
delay_est2 = maxind2 - M; 

%% Ƶ����㻥��أ������ӳ���
fft_x11 = fft(x11,200);
fft_x12 = fft(x12,200);

xc_sp1 = ifft(fft_x11 .* conj(fft_x12));
xc_sp2 = ifft(fft_x12 .* conj(fft_x11));

figure
subplot 211
plot(real(xc_sp1))
axis tight
subplot 212
plot(real(xc_sp2))
axis tight

[~,maxind_freq1] = max(xc_sp1);
delay_est1_freq = maxind_freq1 - 1;

[~,maxind_freq2] = max(xc_sp2);
%% �ӳ�������λӰ���������Ե�
% 236��-20����λӰ����һ�µģ�����Ϊnfft����
% �ӳ�����nfft��Ƚ�Сʱ������ط�ֵ����Խ׼ȷ����ֵԽͻ��
% �ӳ�����nfft�൱���߸���ʱ������ط�ֵ���Ʋ�׼ȷ
delay_est2_freq = maxind_freq2 - 1;
