% 重叠处理后的信号恢复
clear all
close all
clc
%% 读取语音数据
[audio,fs] = audioread('stereo.wav');
[M,N] = size(audio);
frame_time = 10; % 每帧时长
frame_len = fs * frame_time / 1000; % 每帧长度
frame_shift = 128; % 连续帧之间的帧移
lapped_len = frame_len - frame_shift; % 连续帧间重叠点数
frame_num = floor( ( M - lapped_len ) / frame_shift ); % 帧数
mono_audio = audio(:,1);
nfft = 2 ^ nextpow2( frame_len );
freq_num = 1 + nfft / 2;
wint = hanning(frame_len);
out = [];
for i = 1 : frame_num
    iniind = 1 + ( i - 1 ) * frame_shift;
    endind = iniind + frame_len - 1;
    xt = mono_audio( iniind : endind ); % 分帧
    xt_win = wint .* xt; % 加窗
    fft_xt_win = fft( xt_win, nfft ); % FFT
    pf = fft_xt_win( 1 : freq_num ); % 取半轴频率
    % 补全频谱
    full_pf = [ pf; conj( pf( end - 1 : - 1 : 2 ) ) ];
    % 傅里叶反变换回时域
    inv_full_pf = ifft(full_pf, nfft);
    inv_full_pf_dewin = inv_full_pf( 1 : frame_len ) ./ wint; % 取有效数据段
    if ( i < frame_num )
        out_xt = inv_full_pf_dewin( 1 : frame_shift );
    else
        out_xt = inv_full_pf_dewin;
    end
    out = [ out; out_xt ];
end
Mo = length( out );
t = ( 1 : M ) / fs;
to = ( 1 : Mo ) / fs;
figure
subplot 211
hold on 
plot(t,mono_audio,'b.');
plot(to,out,'r')
xlim([min(t) max(t)])
legend('原始信号','加窗重叠处理的恢复信号')

subplot 212
hold on 
plot(t,mono_audio,'b.');
plot(to,out,'r')
xlim([4.1 4.15])
legend('原始信号','加窗重叠处理的恢复信号')

