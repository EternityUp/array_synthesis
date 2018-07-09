clear all
close all
clc
%% ɢ�����������ϵ��
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
title('ɢ����������ɺ�����Ƶ�ʵĹ�ϵ')


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
title('ɢ����������ɺ��������Ĺ�ϵ')

fd = 0 : 1 : 1000;
cf = sinc( 2 * pi * fd / c );
figure
plot(fd,cf)
grid on
axis tight
xlabel('fd/(hz*m)')
ylabel('Real Coherence')
title('ɢ����������ɺ�����(Ƶ��*����)�Ĺ�ϵ')