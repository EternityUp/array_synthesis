clear all
close all
clc
load vad_decision.mat
load post_snr_course.mat
frame_num = size(vad_decision);
%% 读取语音数据
[audio,fs] = audioread('stereo.wav');
[M,N] = size(audio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_time = 10; % 每帧时长
frame_len = fs * frame_time / 1000; % 每帧长度
frame_shift = 128; % 连续帧之间的帧移
lapped_len = frame_len - frame_shift; % 连续帧间重叠点数
nfft = 2 * frame_len - 1; % 傅里叶变换点数
freq_num = frame_len; % 半轴频率点数
wint = hanning( frame_len ); % 窗函数
c_speed = 343; % 声速
mic_d = 0.05; % 阵元间距
f = ( 0 : nfft / 2 ) * fs / nfft;
t = ( 1 : M ) / fs;
frame_t = ( 0 : ( frame_num - 1 ) ) * frame_shift / fs;

for i = 1 : frame_num - 1
    iniind = 1 + ( i - 1 ) * frame_shift;
    endind = iniind + frame_len - 1;
    if ( endind <= M )
        xin = audio( iniind : endind, : ); % 分帧
    else
        xin = audio( iniind : end, : ); % 分帧
    end
    xin_win = xin .* repmat( wint, 1, N );
    fft_xin_win = fft(xin_win, nfft, 1);
    snr = squeeze(post_snr_course(i,:,:));
    full_snr = [ snr ; snr( end : -1 : 2, : ) ]; % 补齐信噪比
    vad_framei = vad_decision(i,:);
    
    if ( all( vad_framei ) )
        %% gcc-phat
        S12 = fft_xin_win(:,1) .* conj( fft_xin_win(:,2) );
        C12 = fftshift(ifft(full_snr(:,1) .* full_snr(:,2) .* S12 ./ abs(S12)));
        [maxval1,max_ind1] = max(C12);
        delay1 = max_ind1 - frame_len;
        doa1 = -asin( delay1 * c_speed / fs / mic_d );
        %% xcorr
        [C12_,lags] = xcorr(xin_win(:,1),xin_win(:,2));
        [maxval2,max_ind2] = max(C12_);
        delay2 = max_ind2 - frame_len;
        doa2 = -asin( delay2 * c_speed / fs /mic_d );
    else
        doa1 = nan;
        doa2 = nan;
    end
    doa1_course(i) = doa1;
    doa2_course(i) = doa2;
end

figure
hold on
plot(C12 / max(C12))
plot(C12_ / max(C12_),'r')
axis tight



