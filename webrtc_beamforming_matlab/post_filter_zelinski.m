clear all
close all
clc

[audio, fs] = audioread('bf_out_nolap.wav');
[M,N] = size(audio);
%% 音频流读取基本参数
chunk_length_ = 256;
frame_shift = 128;
lapped_len = chunk_length_ - frame_shift;
wint = hamming(chunk_length_);
Hwin = repmat(wint,1,2);
%% 音频长度和帧数
L = length(audio);
frame_num = floor( ( L - lapped_len ) / frame_shift );
freq = ( 0 : chunk_length_ / 2 ) * fs / chunk_length_;
freqbins = length(freq);

out_buf_course = zeros(( frame_num - 1 ) * frame_shift...
    + chunk_length_,1);
Nc = N * ( N - 1 ) / 2;
Cp_Prev = zeros(freqbins,Nc);
Ap_chs_Prev = zeros(freqbins,N);
alpha = 0.8;
tic
for i = 1 : frame_num
    ini_ind = 1 + ( i - 1 ) * frame_shift;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind, : );
    
    
    %% FFT变换
    fft_in_block = fft(buf.*Hwin);
    fft_in_block = fft_in_block(1:freqbins,:);
    %% 后置滤波
    % (信号+噪声)自功率谱估计
    Ap_chs = alpha * Ap_chs_Prev + ...
        ( 1 - alpha ) * fft_in_block .* conj(fft_in_block);
    Ap_chs_Prev = Ap_chs;
    Ap = mean(Ap_chs,2);
    
    % zelinski postfilter
    k = 1;
    for k1 = 1 : N - 1
        for k2 = k1 + 1 : N
            Cp(:,k) = real(fft_in_block(:,k1) .* conj(fft_in_block(:,k2)));
            k = k + 1;
        end
    end
    Cp = alpha * Cp_Prev + ( 1 - alpha ) * Cp;
    Cp_Prev = Cp;
    avgCp = mean(Cp,2);
    
    Hz = avgCp ./ Ap;
    Hz = max(Hz,0.05);
    Hz = min(Hz,1);
    fft_out_block = sum(fft_in_block,2);
    fft_out_block = Hz .* fft_out_block;
    
    %% 构造复数全谱
    full_com_sp = [fft_out_block;conj(fft_out_block(end-1:-1:2,:))];
    output_block_ = (ifft(full_com_sp,[],1));
    avg_output_block = mean(output_block_,2);
    
    out_buf_course(audio_ind) = out_buf_course(audio_ind) + ...
        avg_output_block;
    
end
toc
out_buf_course = real(out_buf_course);
norm_factor = max(abs(out_buf_course));
out_buf_course = out_buf_course / norm_factor;

audiowrite('out_pf_test.wav',out_buf_course,fs);


