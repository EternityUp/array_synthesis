clear all
close all
clc

% [audio,fs] = audioread('s6_section1.wav');
[audio,fs] = audioread('stereo.wav');
audio = audio(:,2);
%% 音频流读取基本参数
chunk_length_ = 160;
frame_offset_ = 0;
block_size_ = 256;
freqbins = 129;
frame_shift = 128;
shift_amount_ = 128;
initial_delay_ =  block_size_ - gcd( chunk_length_, shift_amount_ );
element_count = chunk_length_ + initial_delay_;
input_buffer_ = zeros( element_count, 1 );
output_buffer_ = zeros( element_count, 1 );
input_block_ = zeros( block_size_, 1 );
output_block_ = zeros( block_size_, 1 );
read_pos = 160;
write_pos = 0;
rw_wrap = 1; % 0:SAME_WRAP   1:DIFF_WRAP
%% 音频长度和帧数
L = length(audio);
frame_num = floor( L / chunk_length_ );
%% vad相关参数
blocks = 0;
wint = KaiserBesselDerived(1.5,block_size_);
NIS = 5;
omega = 3; % 频域平滑窗口大小,取奇数
pf_smooth_win = hanning( omega );
% 窗函数能量归一化
pf_smooth_win = pf_smooth_win ./ sqrt( pf_smooth_win' * pf_smooth_win );
alpha_pre_emphasis = 0.9375;
alpha_bkg = 0.8;
alpha_pf = 0.9;
frequency_start = 50;
frequency_end = 4000;
alpha_log = 2;
beyond_ratio_energy = 5;
beyond_ratio_ee = 3;
non_speech_to_speech_frame_len = 20;
speech_to_non_speech_frame_len = 20;
non_speech_to_speech_list = zeros(non_speech_to_speech_frame_len,1);
speech_to_non_speech_list = ones(speech_to_non_speech_frame_len,1);
non_speech_to_speech_ret = 0;
speech_to_non_speech_ret = 0;

is_speech_frame = 0;
is_speech_triggered = 0;
speech_start_threshold = 0.2;
speech_finished_threshold = 0.2;

frequency_start_ind = floor(frequency_start * block_size_ / fs)+1;
frequency_end_ind = floor(frequency_end * block_size_ / fs)+1;
tmp_bkg_sp = zeros(freqbins,1);
smoothed_tmp_bkg_sp = zeros(freqbins,1);
InstantSmoothedPf = zeros(freqbins,1);

is_speech_frame_course = [];
is_speech_triggered_course = [];

eeratio_course = [];
energy_course = [];

ref_eeratio_course = [];
ref_energy_course = [];

tic
for i = 1 : frame_num
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind);
    [input_buffer_,read_pos,write_pos,rw_wrap] = ...
        InputBufferWrite(buf,input_buffer_,chunk_length_,...
        read_pos,write_pos,rw_wrap,element_count);
    first_frame_in_block = frame_offset_;
    blks = 0;
    AvgInFdPf = zeros(freqbins,1);
    while ( first_frame_in_block < chunk_length_ )
        [input_block_,read_pos,write_pos,rw_wrap] = ...
            InputBufferRead(input_buffer_,input_block_,block_size_,...
            read_pos,write_pos,rw_wrap,element_count);
        moved_frames = -block_size_ + shift_amount_;
        [read_pos,write_pos,rw_wrap] = ...
            InputBufferMoveReadPositionBackward(moved_frames,read_pos,...
            write_pos,rw_wrap,element_count);
        %% 预加重
        input_block_ = [ input_block_(1); ... 
            input_block_(2:end) - alpha_pre_emphasis * input_block_(1:end-1) ];
        %% 加窗
        input_block_ = input_block_ .* wint;
        fft_in_block = fft(input_block_);
        in_fd_block = abs(fft_in_block(1:freqbins)) .^ 2;
        ang_in_fd_block = angle(fft_in_block(1:freqbins));
        if ( i <= NIS )
            tmp_bkg_sp = tmp_bkg_sp + in_fd_block;
            blocks = blocks + 1;
        else
            if ( i == NIS + 1 )
                tmp_bkg_sp = tmp_bkg_sp / blocks;
                smoothed_tmp_bkg_sp = smoothPf( tmp_bkg_sp, omega,...
                    pf_smooth_win , freqbins, 1);
                bkg_noise_energy = mean(smoothed_tmp_bkg_sp...
                    (frequency_start_ind:frequency_end_ind));
                bkg_noise_energy_entropy_ratio = ...
                    ComputeLogEnergyEntropyRatio(smoothed_tmp_bkg_sp,...
                    frequency_start_ind,frequency_end_ind,alpha_log);
            end 
            AvgInFdPf = AvgInFdPf + in_fd_block;

        end
        first_frame_in_block = first_frame_in_block + shift_amount_;
        blks = blks + 1;
    end
    if( i <= NIS )
        is_speech_frame_course = [is_speech_frame_course is_speech_frame];
        is_speech_triggered_course = [is_speech_triggered_course is_speech_triggered];
    end
    
    
    if ( i >= NIS + 1 )
        %% 递归平滑以减小波动
        AvgInFdPf = AvgInFdPf / blks;       
        InstantSmoothedPf = alpha_pf * InstantSmoothedPf + ( 1 - alpha_pf ) * ...
            smoothPf( AvgInFdPf, omega,pf_smooth_win , freqbins, 1);
%         InstantSmoothedPf = smoothPf( AvgInFdPf, omega,pf_smooth_win , freqbins, 1);
        instant_energy = mean(InstantSmoothedPf...
            (frequency_start_ind:frequency_end_ind));
        instant_energy_entropy_ratio = ComputeLogEnergyEntropyRatio(InstantSmoothedPf,...
            frequency_start_ind,frequency_end_ind,alpha_log);
        
        energy_course = [energy_course instant_energy];
        eeratio_course = [eeratio_course instant_energy_entropy_ratio];
        
        if (instant_energy > beyond_ratio_energy * bkg_noise_energy)
            is_speech_frame = 1;
        else
            is_speech_frame = 0;
        end
        
        
%         if (instant_energy_entropy_ratio > beyond_ratio_ee * bkg_noise_energy_entropy_ratio)
%             is_speech_frame = 1;
%         else
%             is_speech_frame = 0;
%         end
        
        is_speech_frame_course = [is_speech_frame_course is_speech_frame];
        
%         if(~is_speech_frame)                 
%             smoothed_tmp_bkg_sp = alpha_bkg * smoothed_tmp_bkg_sp + ...
%                 (1-alpha_bkg) * InstantSmoothedPf;
%             bkg_noise_energy = mean(smoothed_tmp_bkg_sp...
%                 (frequency_start_ind:frequency_end_ind));
%             bkg_noise_energy_entropy_ratio = ...
%                 ComputeLogEnergyEntropyRatio(smoothed_tmp_bkg_sp,...
%                 frequency_start_ind,frequency_end_ind,alpha_log);
%         end

        
       if(~is_speech_triggered)
           min_ind = min(i,non_speech_to_speech_frame_len);
           non_speech_to_speech_ret = mod(i,non_speech_to_speech_frame_len);
           
           if(non_speech_to_speech_ret == 0)
               non_speech_to_speech_ret = non_speech_to_speech_frame_len;
           end
           
           non_speech_to_speech_list(non_speech_to_speech_ret) = is_speech_frame;
           sum_true = sum(non_speech_to_speech_list(1:min_ind));
           if( sum_true / min_ind > speech_start_threshold)
               i
               is_speech_triggered = 1;
               non_speech_to_speech_list = zeros(non_speech_to_speech_frame_len,1);
           end
       else
            min_ind2 = min(i,speech_to_non_speech_frame_len);
            speech_to_non_speech_ret = mod(i,speech_to_non_speech_frame_len);
            
            if(speech_to_non_speech_ret == 0)
               speech_to_non_speech_ret = speech_to_non_speech_frame_len;
            end
           
            speech_to_non_speech_list(speech_to_non_speech_ret) = is_speech_frame;
            sum_true2 = sum(speech_to_non_speech_list(1:min_ind2));
            if(sum_true2 / min_ind2 < speech_finished_threshold)
                is_speech_triggered = 0;
                speech_to_non_speech_list = ones(speech_to_non_speech_frame_len,1);
            end
       end
       
       is_speech_triggered_course = [is_speech_triggered_course is_speech_triggered];
       
       ref_energy_course = [ref_energy_course bkg_noise_energy];
       ref_eeratio_course = [ref_eeratio_course bkg_noise_energy_entropy_ratio];
    end
    frame_offset_ = first_frame_in_block - chunk_length_;
end
toc



figure
plot(1:frame_num,is_speech_frame_course,'bo')
hold on
plot(1:frame_num,is_speech_triggered_course,'r*')
ylim([-1 2])
xlim([ 1 frame_num])


figure
subplot 211
plot(audio)
axis tight
subplot 212
plot(energy_course)
hold on
plot(ref_energy_course,'r')
plot(beyond_ratio_energy * ref_energy_course,'g')

figure
subplot 211
plot(audio)
axis tight
subplot 212
plot(eeratio_course,'b')
hold on
plot(ref_eeratio_course,'r')
plot(beyond_ratio_ee * ref_eeratio_course,'g')






