clear all
close all
clc
N = 1024;
x1 = randn( N ,1 );
M = 256;
delay = 20;
x11 = x1( 1 : M );
x12 = x1( 1 + delay : M + delay );
%% 时域估计互相关  估计延迟量
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

%% 频域计算互相关，估计延迟量
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
%% 延迟量的相位影响是周期性的
% 236与-20的相位影响是一致的，周期为nfft长度
% 延迟量与nfft相比较小时，互相关峰值估计越准确，峰值越突出
% 延迟量与nfft相当或者更大时，互相关峰值估计不准确
delay_est2_freq = maxind_freq2 - 1;
