clear all
close all
clc
[audio,fs] = audioread('cos_stereo.wav');
[M,N] = size(audio);
%% 音频流读取基本参数
chunk_length_ = 320;
%% 音频长度和帧数
L = length(audio);
frame_num = floor( L / chunk_length_ );
%% 其他参数
theta = ( 0 : 10 : 180 ) * pi / 180;
Lt = length(theta);
freq = ( 0 : chunk_length_ / 2 ) * fs / chunk_length_;
freqbins = length(freq);
%% 阵型分析
array_shape = [0 0 0;0.05 0 0];
centered = mean(array_shape,1);
centered_array = array_shape - repmat(centered,N,1);
% centered_array = array_shape;
doa = 60 * pi / 180;
sound_path = zeros(N,1);
sound_speed = 340;
% 计算声程
for k = 1 : N
    sound_path(k) = centered_array(k,1) * cos(doa) + ...
        centered_array(k,2) * sin(doa);
end
tao = sound_path / sound_speed;

fft_in_block = zeros(freqbins,N);
fft_out_block = zeros(freqbins,N);

out_buf_course = [];

tic
for i = 1 : frame_num
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind, : );
    
    %% FFT变换
    fft_in_block = fft(buf);
    fft_in_block = fft_in_block(1:freqbins,:);
    %% 常规波束形成
    for j = 1 : freqbins
        fm = freq(j);
        sv = exp( -1j * 2 * pi * fm * ( tao - tao(1) ) );
        xfk = fft_in_block(j,:).'; % 单频点谱值向量
        fft_out_block(j,:) = conj(sv) .* xfk;
    end
    %% 构造复数全谱
    full_com_sp = [fft_out_block;conj(fft_out_block(end-1:-1:2,:))];
    output_block_ = real(ifft(full_com_sp,[],1));
    out_buf_course = [out_buf_course;output_block_];
end
toc

norm_factor = max(max(abs(out_buf_course)));
out_buf_course = out_buf_course / norm_factor;

audiowrite('bf_out_nolap.wav',out_buf_course,fs);

