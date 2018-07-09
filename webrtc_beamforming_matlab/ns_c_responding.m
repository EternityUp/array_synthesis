clear all
close all
clc

% [audio,fs] = audioread('s6_section1.wav');
[audio,fs] = audioread('cos_stereo.wav');
audio = audio(:,1);
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
%% 降噪相关参数
blocks = 0;
wint = KaiserBesselDerived(1.5,block_size_);
out_buf_course = [];
NIS = 5;
omega = 3; % 频域平滑窗口大小,取奇数
pf_smooth_win = hanning( omega );
alpha_s = 0.8; % 功率谱平滑因子
NL = 125;
delta = 5; % 噪声估计的似然比阈值
alpha_p = 0.8; % 估计有语音条件概率的平滑因子
alpha_d = 0.95; % 噪声谱估计的平滑因子
beta = 0.7;        % 估计递归平均先验信噪比的平滑因子
% 窗函数能量归一化
pf_smooth_win = pf_smooth_win ./ sqrt( pf_smooth_win' * pf_smooth_win );
omega_local = 3;   % 局部平滑窗口大小
h_local = hanning(omega_local);  % 局部平滑汉宁窗
h_local = h_local ./ sqrt( h_local' * h_local ); % 窗函数能量归一化
omega_global = 31; % 全局平滑窗口大小
h_global = hanning(omega_global); % 全局平滑汉宁窗
h_global = h_global ./ sqrt( h_global' * h_global ); % 窗函数能量归一化
sigma_min = -10;
sigma_max = -5;
sigma_pmin = 0;
sigma_pmax = 10;
q_max = 0.95;      % 无语音先验概率的最大值
alpha = 0.99 * ones(freqbins,1); % 用于平衡噪声消除程度与过渡阶段引入的语音失真的平滑因子
PostSNR_Prev  = 10 .^ ( -25 / 20 ) * ones(freqbins,1); % 初始化后验信噪比:幅度
GainH1_Prev   = 10 .^ ( -25 / 20 ) * ones(freqbins,1); % 初始化有语音时的条件增益:幅度
Gmin = -30;   % 无语音时的最小增益值:dB





tmp_bkg_sp = zeros(freqbins,1);
smoothed_tmp_bkg_sp = zeros(freqbins,1);
tmp_Ef_Prev = zeros(freqbins,1);
min_Ef_Prev = zeros(freqbins,1);
NoiseSpectrumEst_Prev = zeros(freqbins,1);
InstantSmoothedPf = zeros(freqbins,1);
RecursiveEf =  zeros(freqbins,1);
RecursiveEf_Prev =  zeros(freqbins,1);
ratio =  zeros(freqbins,1);
SpeechPresenceProb_Prev = 0.01 * ones(freqbins,1);
RecurAvgPriorSNR =  zeros(freqbins,1);
PriorSNR_Prev = 10 .^ ( -25 / 20 ) * ones(freqbins, 1);
RecurAvgPriorSNR_PeakPrev = 10 ^ ( -25 / 20 );
AvgRecurAvgPriorSNRPrev   = 10 ^ ( -25 / 20 );
RecurAvgPriorSNRPrev = 10 .^ ( -25 / 20 ) * ones(freqbins, 1);
output_block_ = zeros(block_size_,1);

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
    while ( first_frame_in_block < chunk_length_ )
        [input_block_,read_pos,write_pos,rw_wrap] = ...
            InputBufferRead(input_buffer_,input_block_,block_size_,...
            read_pos,write_pos,rw_wrap,element_count);
        moved_frames = -block_size_ + shift_amount_;
        [read_pos,write_pos,rw_wrap] = ...
            InputBufferMoveReadPositionBackward(moved_frames,read_pos,...
            write_pos,rw_wrap,element_count);
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
                RecursiveEf_Prev = smoothed_tmp_bkg_sp;
                tmp_Ef_Prev = smoothed_tmp_bkg_sp;
                min_Ef_Prev = smoothed_tmp_bkg_sp;
                NoiseSpectrumEst_Prev = tmp_bkg_sp;
            end
            
            InstantSmoothedPf = smoothPf( in_fd_block, omega,...
                pf_smooth_win , freqbins, 1);
            RecursiveEf = alpha_s * RecursiveEf_Prev + ...
                ( 1 - alpha_s ) * InstantSmoothedPf;
            min_Ef = min( min_Ef_Prev, RecursiveEf );
            tmp_Ef = min( tmp_Ef_Prev, RecursiveEf );
            if ( mod( i, NL ) == 0 )
                min_Ef = min( tmp_Ef_Prev, RecursiveEf );
                tmp_Ef = RecursiveEf;
            end
            ratio = RecursiveEf ./ min_Ef;
            Ikl = ratio > delta;
            SpeechPresenceProb = alpha_p .* SpeechPresenceProb_Prev + ...
                ( 1 - alpha_p ) * Ikl;
            alpha_df = alpha_d + ( 1 - alpha_d ) .* SpeechPresenceProb;
            NoiseSpectrumEst = alpha_df .* NoiseSpectrumEst_Prev + ...
                ( 1 - alpha_df ) .* InstantSmoothedPf;
            RecursiveEf_Prev        = RecursiveEf;
            min_Ef_Prev            = min_Ef;
            tmp_Ef_Prev            = tmp_Ef;
            SpeechPresenceProb_Prev = SpeechPresenceProb;
            NoiseSpectrumEst_Prev   = NoiseSpectrumEst;
            RecurAvgPriorSNR = beta .* RecurAvgPriorSNRPrev + ...
                ( 1 - beta ) .* PriorSNR_Prev;
            % local
            RecurAvgPriorSNR_Local = ...
                smoothPf( RecurAvgPriorSNR, omega_local, h_local , freqbins, 1);
            P_Local = RecurAvgPriorSNR_Local;
            % global
            RecurAvgPriorSNR_Global = ...
                smoothPf( RecurAvgPriorSNR, omega_global, h_global , freqbins, 1);
            P_Global = RecurAvgPriorSNR_Global;
            
            for k = 1 : freqbins
                if( P_Local( k ) <= 10 ^ ( sigma_min / 20 ) )
                    P_Local( k ) = 0;
                elseif( P_Local( k ) >= 10 ^ ( sigma_max / 20 ) )
                    P_Local( k ) = 1;
                else
                    P_Local( k ) = log10( P_Local( k ) / 10 ^ ( sigma_min / 20 ) ) ...
                        / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
                end
                
                if( P_Global( k ) <= 10 ^ ( sigma_min / 20 ) )
                    P_Global( k ) = 0;
                elseif( P_Global( k ) >= 10 ^ ( sigma_max / 20 ) )
                    P_Global( k ) = 1;
                else
                    P_Global( k ) = log10( P_Global( k ) / 10 ^ ( sigma_min / 20 ) ) ...
                        / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
                end
            end
            AvgRecurAvgPriorSNR = mean( RecurAvgPriorSNR( 1 : ( frame_shift / 2 + 1), : ), 1 );
            th1 = 10 ^ ( sigma_min / 20 ) * 10 .^ ( RecurAvgPriorSNR_PeakPrev / 20 );
            th2 = 10 ^ ( sigma_max / 20 ) * 10 .^ ( RecurAvgPriorSNR_PeakPrev / 20 );
            % 计算mu_l
            if( AvgRecurAvgPriorSNR <= th1 )
                mu_l = 0;
            elseif( AvgRecurAvgPriorSNR >= th2 )
                mu_l = 1;
            else
                mu_l = log10( AvgRecurAvgPriorSNR / th1 ) ...
                    / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
            end
            % 计算P_Frame
            if( AvgRecurAvgPriorSNR <= 10 ^ ( sigma_min / 20 ) )
                P_Frame = 0;
            elseif( AvgRecurAvgPriorSNR <= AvgRecurAvgPriorSNRPrev )
                P_Frame = mu_l;
            else
                RecurAvgPriorSNR_Peak = min( max( AvgRecurAvgPriorSNR, 10 ^ ( sigma_pmin / 20 ) ),...
                    10 ^ ( sigma_pmax / 20 ) );
                RecurAvgPriorSNR_PeakPrev = RecurAvgPriorSNR_Peak;
                P_Frame = 1;
            end
            % 计算无语音概率
            SpeechAbsenceProb = 1 - ...
                P_Local .* P_Global .* P_Frame;
            % 限制无语音概率最大值为0.95
            SpeechAbsenceProb( SpeechAbsenceProb >= q_max ) = q_max;
            % ------------更新相关变量-------------------%
            AvgRecurAvgPriorSNRPrev = AvgRecurAvgPriorSNR;
            RecurAvgPriorSNRPrev = RecurAvgPriorSNR;
            % 计算后验信噪比
            PostSNRCurr = in_fd_block ./ NoiseSpectrumEst; % 幅度
            
            PriorSNRCurr = alpha .* ( GainH1_Prev .^ 2 ) .* PostSNR_Prev + ...
                ( 1 - alpha ) .* max( PostSNRCurr - 1, 0 );
            vkPost_Prior = PriorSNRCurr .* PostSNRCurr ./ ( 1 + PriorSNRCurr );
            x = expint( vkPost_Prior );
            GainH1_Curr = PriorSNRCurr ./ ( 1 + PriorSNRCurr ) .* ...
                exp( expint( vkPost_Prior ) / 2 );
            PriorSNR_Prev = PriorSNRCurr;
            PostSNR_Prev  = PostSNRCurr;
            GainH1_Prev   = GainH1_Curr;
            GainH1_Prev   = max(GainH1_Curr,10 .^ ( Gmin / 20 ));
            pkl = ( 1 + SpeechAbsenceProb ./ ( 1 - SpeechAbsenceProb ) ...
                .* ( 1 + PriorSNRCurr ) .* exp( -vkPost_Prior ) );
            pkl = 1 ./ pkl;
            Gainf = ( GainH1_Curr .^ pkl ) .* ...
                ( 10 ^ ( Gmin / 20 ) .^ ( 1 - pkl ) );
            out_fd_block = Gainf .^ 2 .* in_fd_block;
            % 构造复数谱
            half_com_sp = sqrt(out_fd_block) .* exp(1j*ang_in_fd_block);
            full_com_sp = [half_com_sp;conj(half_com_sp(end-1:-1:2))];
            output_block_ = real(ifft(full_com_sp));
        end
        output_block_ = output_block_ .* wint;
        % AddFrames
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
    end
    
    % Copy output buffer to output
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
    frame_offset_ = first_frame_in_block - chunk_length_;
    out_buf_course = [out_buf_course;out_buf];
end
toc
out_buf_course = out_buf_course / max(abs(out_buf_course));

figure
plot(audio)
hold on
plot(out_buf_course(225:end),'r')
grid on
axis tight
legend('模拟输入序列','未经算法处理的输出序列')



audiowrite('signle_ns.wav',out_buf_course,fs);

















