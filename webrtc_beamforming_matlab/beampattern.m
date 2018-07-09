clear all
close all
clc

M = 4; % ��Ԫ��Ŀ
f = 6000; % Ƶ��
c = 340; % ����
lambda = c / f; %����
d = 0.05; % ������Ԫ���

theta0 = 0 * pi / 180; % �����
theta = ( -90 : 0.1 : 90 ) * pi / 180; % ɨ���
Lt = length(theta);

pp = zeros(Lt,1);
mu = pi * f * d / c;
for i = 1 : Lt
    sindf = sin(theta(i)) - sin(theta0);
    if(sindf == 0)
        pp(i) = M ^ 2;
    else
        pp(i) = sin( M * mu * sindf ) ^ 2 / sin( mu * sindf ) ^ 2;
    end
end

ppdB = 10 * log10(pp);
ppdB = ppdB - max(ppdB);

figure
plot(theta*180/pi,ppdB)
xlabel('ɨ��Ƕ�/��')
ylabel('Am/dB')
grid on
axis tight






