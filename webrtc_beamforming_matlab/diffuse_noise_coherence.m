clear all
close all
clc
%% 散射噪声场相关系数
c = 340;
f = 0 : 10 : 8000;
dij = 0.1;
cf = sinc( 2 * pi * f * dij / c );
figure
plot(f,cf)
grid on
axis tight
xlabel('f/Hz')
ylabel('Real Coherence')
title('散射噪声场相干函数与频率的关系')


f = 0 : 10 : 8000;
dij = 0.01:0.01:1;
cf = sinc( 2 * pi * f' * dij / c );
figure
imagesc(f,dij,cf');
xlabel('f/Hz')
ylabel('distance/m')
figure
plot(dij,cf(100,:))
grid on
axis tight
xlabel('distance/m')
ylabel('Real Coherence')
title('散射噪声场相干函数与距离的关系')

fd = 0 : 1 : 1000;
cf = sinc( 2 * pi * fd / c );
figure
plot(fd,cf)
grid on
axis tight
xlabel('fd/(hz*m)')
ylabel('Real Coherence')
title('散射噪声场相干函数与(频率*距离)的关系')