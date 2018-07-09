clear all
close all
clc
% load Array_Xt.mat
[audio,fs] = audioread('cos_stereo.wav');
% audio = Array_Xt';
% audio = audio(2.71e4:1:3e4,:);
[M,N] = size(audio);
%% 音频流读取基本参数
chunk_length_ = 160;
frame_offset_ = 0;
block_size_ = 256;
freqbins = 129;
frame_shift = 128;
shift_amount_ = 128;
initial_delay_ =  block_size_ - gcd( chunk_length_, shift_amount_ );
element_count = chunk_length_ + initial_delay_;
input_buffer_ = zeros( element_count, N );
output_buffer_ = zeros( element_count, N );
input_block_ = zeros( block_size_, N );
output_block_ = zeros( block_size_, N );
read_pos = 160;
write_pos = 0;
rw_wrap = 1; % 0:SAME_WRAP   1:DIFF_WRAP
%% 音频长度和帧数
L = length(audio);
frame_num = floor( L / chunk_length_ );
%% 其他参数
wint = KaiserBesselDerived(1.5,block_size_);
frequency_start = 50;
frequency_end = 4000;
frequency_start_ind = floor(frequency_start * block_size_ / fs) + 1;
frequency_end_ind = floor(frequency_end * block_size_ / fs) + 1;
theta = ( 0 : 10 : 180 ) * pi / 180;
Lt = length(theta);
alpha_smooth = 0.8;
freq = ( 0 : block_size_ / 2 ) * fs / block_size_;

%% 阵型分析
% array_shape = [0 0 0;0.05 0 0;0.1 0 0;0.15 0 0; ...
%     0.20 0 0;0.25 0 0;0.30 0 0;0.35 0 0];
% d = ( 0 : N - 1 ) * 0.05;
array_shape = [0 0 0;0.05 0 0];
centered = mean(array_shape,1);
centered_array = array_shape - repmat(centered,N,1);
sound_path = zeros(N,1);
sound_speed = 340;

tic
doa = zeros(frame_num,1);
I = ones(N,1);
for i = 1 : frame_num
    i
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind, : );
    [input_buffer_,read_pos,write_pos,rw_wrap] = ...
            InputBufferWrite(buf,input_buffer_,chunk_length_,...
            read_pos,write_pos,rw_wrap,element_count);
    first_frame_in_block = frame_offset_;
    blks = 0;
    fft_in_block= zeros(freqbins,N);
    theta_est = zeros(2,1);
    while ( first_frame_in_block < chunk_length_ )
        Pf_out = zeros(Lt,1);
        [input_block_,read_pos,write_pos,rw_wrap] = ...
            InputBufferRead(input_buffer_,input_block_,block_size_,...
            read_pos,write_pos,rw_wrap,element_count);
        [read_pos,write_pos,rw_wrap] = ...
            InputBufferMoveReadPositionBackward(-shift_amount_,read_pos,...
            write_pos,rw_wrap,element_count);
        %% 加窗
        input_block_ = input_block_ .* repmat(wint,1,N);
        
        %% FFT变换
        fft_in_block = fft(input_block_);
        fft_in_block = fft_in_block(1:freqbins,:);
        in_fd_block = abs(fft_in_block) .^ 2;
        ang_in_fd_block = angle(fft_in_block);
        
        %% 波达方向估计
        
        for j = 1 : Lt
            RSCTM = zeros(N,N);
            Y_stmv = zeros(N,1);
            theta_scan = theta(j);
            % 计算声程
            for k = 1 : N
                sound_path(k) = centered_array(k,1) * cos(theta_scan) + ...
                    centered_array(k,2) * sin(theta_scan);
            end
            tao = sound_path / sound_speed; 
            % 计算阵列流形
            for m = frequency_start_ind : frequency_end_ind;
                fm = freq(m);
                sv = exp( -1j * 2 * pi * fm * tao );
                xfk = fft_in_block(m,:).'; % 单频点谱值向量
                %% cbf
                Pf_out(j) = Pf_out(j) + abs(sv' * ( xfk * xfk' ) * sv);
                %% mvdr
%                 Rxf_mvdr = xfk * xfk' + 0.01 * eye( N );
%                 Pf_out(j) = Pf_out(j) + 1 / abs( sv' / Rxf_mvdr * sv );
                %% stmv
%                 Tf = diag(sv);
%                 Rxf_stmv = xfk * xfk';
%                 RSCTM = RSCTM + Tf' * Rxf_stmv * Tf;
            end
%             RSCTM = RSCTM + 0.01 * eye( N );
%             wo = ( RSCTM \ I ) / ( I' / RSCTM * I );
%             Pf_out(j) = 1 / abs( I' / RSCTM * I );

        end
        [max_pf,max_index] = max(Pf_out);
        first_frame_in_block = first_frame_in_block + shift_amount_;
        blks = blks + 1;
        theta_est(blks) = theta(max_index) * 180 / pi;
    end
    doa(i) = sum(theta_est(1:blks)) / blks;
    frame_offset_ = first_frame_in_block - chunk_length_;
end
toc

figure
plot(freq,in_fd_block);
xlabel('f/Hz')
ylabel('Am')
grid on
axis tight


Pf_out_dB = 10*log10(Pf_out);
Pf_out_dB = Pf_out_dB - max(Pf_out_dB);
figure
plot(theta*180/pi,Pf_out_dB);
grid on
axis tight

figure
plot(1:frame_num,doa,'*-')
xlabel('帧次')
ylabel('方位估计值')
xlim([1 frame_num])
grid on
