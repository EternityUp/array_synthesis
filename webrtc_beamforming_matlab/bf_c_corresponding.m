clear all
close all
clc
[audio,fs] = audioread('cos_stereo.wav');
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
theta = ( 0 : 10 : 180 ) * pi / 180;
Lt = length(theta);
alpha_smooth = 0.8;
freq = ( 0 : block_size_ / 2 ) * fs / block_size_;
%% 阵型分析
array_shape = [0 0 0;0.05 0 0];
centered = mean(array_shape,1);
centered_array = array_shape - repmat(centered,N,1);
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
        
        %% 加窗
        output_block_ = output_block_ .* repmat(wint,1,2);
        %% AddFrames
        first_frame_in_blocki = first_frame_in_block;
        result_ind_range = (first_frame_in_blocki + 1) : ...
            (first_frame_in_blocki + block_size_);
        b_ind_range = 1 : block_size_;
        a_ind_range = (first_frame_in_blocki + 1) : ...
            (first_frame_in_blocki + block_size_);
        output_buffer_(result_ind_range, :) =  ...
            output_buffer_(a_ind_range, :) + ...
            output_block_(b_ind_range,:);
        
        
        first_frame_in_block = first_frame_in_block + shift_amount_;
        blks = blks + 1;
    end
    
    %% Copy output buffer to output
    % CopyFrames
    out_buf = output_buffer_(1:chunk_length_,:);
    % MoveFrames ----> output_buffer
    dst_ind_range = 1 : initial_delay_;
    src_ind_range = ( chunk_length_ + 1 ) : ...
        ( chunk_length_ + initial_delay_ );
    output_buffer_(dst_ind_range,:) = output_buffer_(src_ind_range,:);
    % ZeroOut ----> output_buffer
    src_ind_range = ( initial_delay_ + 1 ) : ...
        ( chunk_length_ + initial_delay_ );
    output_buffer_(src_ind_range,:) = 0;
    
    out_buf_course = [out_buf_course;out_buf];
    
    frame_offset_ = first_frame_in_block - chunk_length_;
end
toc

norm_factor = max(max(abs(out_buf_course)));
out_buf_course = out_buf_course / norm_factor;

audiowrite('bf_out.wav',out_buf_course,fs);

