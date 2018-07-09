clear all
close all
clc
M = 8;
N = 256;
s = 0.5 * sin( 2 * pi * 2 * ( 0 : N - 1 ) / 30 );
noise = randn(M,N);
x = repmat(s,M,1) + noise;
Nc = M * ( M - 1 ) / 2;
%% 计算各通道噪声频谱以及通道间噪声相干函数
fft_n = fft(noise,[],2);
pfn = fft_n .* conj(fft_n);
avg_pfn = mean(pfn); % 非相干叠加
cfn = zeros(Nc,N);
k = 1;
for i = 1 : M - 1
    for j = i + 1 : M
        cfn(k,:) = fft_n(i,:) .* conj(fft_n(j,:)) ./ avg_pfn;
        % 下面计算方式不正确，因为没有考虑相干谱估计的统计特性
%         cfn(k,:) = fft_n(i,:) .* conj(fft_n(j,:)) ./ ...
%             sqrt((fft_n(i,:) .* conj(fft_n(i,:))) .* ...
%             (fft_n(j,:) .* conj(fft_n(j,:))));
        k = k + 1;
    end
end
norm_factor = max(max(abs(cfn)));
cfn = cfn / norm_factor;

%% 计算各通道频谱,估计自相关功率谱
fft_x = fft(x,[],2);
Ap = mean(fft_x .* conj(fft_x),1);

%% 估计互相关功率谱
Cp = zeros(Nc, N);
kk = 1;
for i = 1 : M - 1
    for j = i + 1 : M
        xc_ij = real(fft_x(i,:) .* conj(fft_x(j,:)));
        xa_i = fft_x(i,:) .* conj(fft_x(i,:));
        xa_j = fft_x(j,:) .* conj(fft_x(j,:));
        cfn_real_ij = real(cfn(kk,:));
        cfn_real_ij( cfn_real_ij >= 1 ) = 0.95; % 相干系数阈值
        Cp(kk,:) = ( xc_ij - 0.5 * cfn_real_ij .* ( xa_i + xa_j ) ) ...
            ./ ( 1 - cfn_real_ij );
        kk = kk + 1;
    end
end
avgCp = mean(Cp,1);

Hz = avgCp ./ Ap;

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
title('single channel mccowan output')

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
title('multiple channels mccowan output')


figure
plot(out_x,'b-o')
hold on
plot(out_x2,'r-*')
axis tight
legend('single channel mccowan output','multiple channels mccowan output')


