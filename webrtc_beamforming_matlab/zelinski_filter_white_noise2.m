close all
clear all
clc
M = 8;
N = 256;
s = 0.5 * sin( 2 * pi * 2 * ( 0 : N - 1 ) / 30 );
x = repmat(s,M,1) + randn(M,N);
load x.mat
%% 估计自相关功率谱
fft_x = fft(x,[],2);
Ap = mean(fft_x .* conj(fft_x),1);
%% 估计互功率谱
Nc = M * ( M - 1 ) / 2;
Cp = zeros(Nc, N);
Sk = zeros(Nc, N);
Vk = zeros(Nc, N);
k = 1;
for i = 1 : ( M - 1 )
    for j = ( i + 1 ) : M
        Cp(k,:) = real( fft_x(i,:) .* conj(fft_x(j,:)) );
        Sk(k,:) = Cp(k,:);
        Sk(k,(Sk(k,:) < 0)) = 0;
        Vk(k,:) = Cp(k,:);
        Vk(k,(Vk(k,:) > 0)) = 0;
        k = k + 1;
    end
end

b_Sk = Sk > 0;
b_Sk_Num = sum(b_Sk,1);
b_Sk_Num(b_Sk_Num == 0) = 1;
b_Vk = Vk < 0;
b_Vk_Num = sum(b_Vk,1);
b_Vk_Num(b_Vk_Num == 0) = 1;

SkEst = sum(Sk,1) ./ b_Sk_Num;

VkEst_Square = sum(Vk.^2,1) ./ b_Vk_Num;

SmoothSkEst = smooth(SkEst,5,'moving');
SmoothVkEst = smooth(VkEst_Square,5,'moving');

alpha_k = SmoothSkEst.^2 ./ ( SmoothSkEst.^2 + SmoothVkEst / Nc );
avgCp = mean(Cp,1);

Cpp = alpha_k' .* avgCp;

Hz = Cpp ./ Ap;

%% 单通道后置滤波
avgx = mean(x,1);
fft_avgx = fft(avgx);
out_x = real( ifft( Hz .* fft_avgx ) );
figure
subplot 311
plot(x(1,:))
axis tight
title('ch1')
subplot 312
plot(avgx)
axis tight
title('beamformer output')
subplot 313
plot(out_x)
axis tight
title('single channel zelinski output')

%% 多通道后置滤波
Hz_chs = repmat(Hz,M,1);
out_x_chs = real( ifft( Hz_chs .* fft_x, [], 2 ) );
out_x2 = mean(out_x_chs,1);
figure
subplot 311
plot(x(1,:))
axis tight
title('ch1')
subplot 312
plot(avgx)
axis tight
title('beamformer output')
subplot 313
plot(out_x2)
axis tight
title('multiple channels zelinski output')


figure
plot(out_x,'b-o')
hold on
plot(out_x2,'r-*')
axis tight
legend('single channel zelinski output','multiple channels zelinski output')


