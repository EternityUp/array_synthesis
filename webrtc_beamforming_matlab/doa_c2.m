clear all
close all
clc
clear all
close all
clc
% load Array_Xt.mat
[audio,fs] = audioread('cos_stereo_seg.wav');
% audio = Array_Xt';
% audio = audio(2.71e4:1:3e4,:);
[M,N] = size(audio);
%% 音频流读取基本参数
chunk_length_ = 160;
%% 音频长度和帧数
L = length(audio);
frame_num = floor( L / chunk_length_ );
%% 其他参数
frequency_start = 300;
frequency_end = 6000;
frequency_start_ind = floor(frequency_start * chunk_length_ / fs) + 1;
frequency_end_ind = floor(frequency_end * chunk_length_ / fs) + 1;
theta = ( 0 : 1 : 180 ) * pi / 180;
Lt = length(theta);
alpha_smooth = 0.8;
freq = ( 0 : chunk_length_ / 2 ) * fs / chunk_length_;
freqbins = chunk_length_ / 2 + 1;
%% 阵型分析
% array_shape = [0 0 0;0.05 0 0;0.1 0 0;0.15 0 0; ...
%     0.20 0 0;0.25 0 0;0.30 0 0;0.35 0 0];
array_shape = [0 0 0;0.05 0 0];
d = ( 0 : N - 1 ) * 0.05;
centered = mean(array_shape,1);
centered_array = array_shape - repmat(centered,N,1);
sound_path = zeros(N,1);
sound_speed = 340;

tic
for i = 5 : 5
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind, : );
    fft_in_block = fft(buf,chunk_length_,1);
    fft_in_block = fft_in_block(1:freqbins,:);
    in_fd_block = abs(fft_in_block) .^ 2;
    ang_in_fd_block = angle(fft_in_block);
    
    theta_pf = zeros(Lt,1);
    %% 波达方向估计
    for j = 1 : Lt
        theta_scan = theta(j);
        % 计算声程
        for k = 1 : N
            sound_path(k) = centered_array(k,1) * cos(theta_scan) + ...
                centered_array(k,2) * sin(theta_scan);
        end
        tao = sound_path / sound_speed;
%         tao = d' * cos(theta_scan) / sound_speed;
        % 计算阵列流形
        for m = frequency_start_ind : frequency_end_ind;
            fm = freq(m);
            sv = exp( -1j * 2 * pi * fm * tao );
            xfk = fft_in_block(m,:).'; % 单频点谱值向量
%             xfk = xfk ./ abs(xfk);
            theta_pf(j) = theta_pf(j) + abs(sv' * ( xfk * xfk' ) * sv);
        end
    end
    
    
    
end
toc

figure
plot(freq,10*log10(in_fd_block));
xlabel('f/Hz')
ylabel('Am')
grid on
axis tight

theta_pf_dB = 10*log10(theta_pf);
theta_pf_dB = theta_pf_dB - max(theta_pf_dB);
figure
plot(theta*180/pi,theta_pf_dB);
grid on
axis tight


%% 互相关求延时
% x1 = buf(:,1);
% x2 = buf(:,2);
% [xc,lags] = xcorr(x1,x2,'unbiased');
% figure
% plot(lags,xc)



