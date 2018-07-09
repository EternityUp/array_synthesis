clear all
close all
clc

M = 4; % 阵元数目
f = 6000; % 频率
c = 340; % 声速
lambda = c / f; %波长
d = 0.05; % 线阵阵元间距

theta0 = 0 * pi / 180; % 入射角
theta = ( -90 : 0.1 : 90 ) * pi / 180; % 扫描角
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
xlabel('扫描角度/°')
ylabel('Am/dB')
grid on
axis tight






