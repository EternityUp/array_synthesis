clear all
close all
clc
M = 2; % 阵元数目
c = 340; % 声速
d = 0.05; % 线阵阵元间距

theta0 = 0 * pi / 180; % 入射角
theta = ( -90 : 0.1 : 90 ) * pi / 180; % 扫描角
Lt = length(theta);

f = 50:50:4000; % 波束形成频率范围
Lf = length(f);
taoO = ( 0 : ( M - 1 ) )' * d * sin( theta0 ) / c;
Pf_cbf = zeros(Lt,1);
%% cbf
tic
for i = 1 : Lt
    thetai = theta(i);
    for j = 1 : Lf
        fj = f(j);
        tao = ( 0 : ( M - 1 ) )' * d * sin( thetai ) / c; 
        sv = exp( -1j * 2 * pi * fj * tao );
        svo = exp ( -1j * 2 * pi * fj * taoO );
        Pf_cbf(i) = Pf_cbf(i) + abs( sv' * svo );  
    end
end
toc

Pf_cbf_dB = 20 * log10(Pf_cbf);
Pf_cbf_dB = Pf_cbf_dB - max(Pf_cbf_dB);
figure
plot(theta*180/pi,Pf_cbf_dB);
xlabel('扫描角度/°')
ylabel('Am/dB')
grid on
axis tight
title('CBF')

%% mvdr
tic
Pf_mvdr = zeros(Lt,1);
for i = 1 : Lt
    thetai = theta(i);
    for j = 1 : Lf
        fj = f(j);
        tao = ( 0 : ( M - 1 ) )' * d * sin( thetai ) / c; 
        sv = exp( -1j * 2 * pi * fj * tao );
        svo = exp ( -1j * 2 * pi * fj * taoO );
        Rxf = svo * svo' + 0.01 * eye( M );
        Pf_mvdr(i) = Pf_mvdr(i) + 1 / abs( sv' / Rxf * sv );
    end
end
toc

Pf_mvdr_dB = 20 * log10(Pf_mvdr);
Pf_mvdr_dB = Pf_mvdr_dB - max(Pf_mvdr_dB);
figure
plot(theta*180/pi,Pf_mvdr_dB);
xlabel('扫描角度/°')
ylabel('Am/dB')
grid on
axis tight
title('MVDR')



%% stmv
tic
Pf_stmv = zeros(Lt,1);
I = ones(M,1);
for i = 1 : Lt
    thetai = theta(i);
    RSCTM = zeros(M,M);
    for j = 1 : Lf
        fj = f(j);
        tao = ( 0 : ( M - 1 ) )' * d * sin( thetai ) / c; 
        sv = exp( -1j * 2 * pi * fj * tao );
        svo = exp ( -1j * 2 * pi * fj * taoO );
        Tf = diag(sv);
        Rxf = svo * svo';
        RSCTM = RSCTM + Tf' * Rxf * Tf;
    end
    RSCTM = RSCTM + 0.01 * eye( M );
    Pf_stmv(i) = 1 / abs( I' / RSCTM * I ); 
end
toc

Pf_stmv_dB = 20 * log10(Pf_stmv);
Pf_stmv_dB = Pf_stmv_dB - max(Pf_stmv_dB);
figure
plot(theta*180/pi,Pf_stmv_dB);
xlabel('扫描角度/°')
ylabel('Am/dB')
grid on
axis tight
title('STMV')
