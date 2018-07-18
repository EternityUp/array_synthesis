clear all
close all
clc

% [audio,fs] = audioread('cos_stereo.wav');
% [audio,fs] = audioread('stereo_linear.wav');
% [audio,fs] = audioread('stereo_circular.wav');
[audio,fs] = audioread('xmos_record1_chs4.wav');
N = size(audio,2);
% array_shape = construct_array_loc(N,0.05,1); % 构造阵列阵元位置
load  array_loc_xmos.mat
array_shape = array_loc_xmos;
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
%% 降噪相关参数
blocks = 0;
wint = KaiserBesselDerived(1.5,block_size_);
NIS = 10;
omega = 3; % 频域平滑窗口大小,取奇数
pf_smooth_win = hanning( omega );
alpha_s = 0.8; % 功率谱平滑因子
NL = 125;
delta = 1; % 噪声估计的似然比阈值
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
alpha = 0.99 * ones(freqbins,N); % 用于平衡噪声消除程度与过渡阶段引入的语音失真的平滑因子
PostSNR_Prev  = 10 .^ ( -25 / 20 ) * ones(freqbins,N); % 初始化后验信噪比:幅度
GainH1_Prev   = 10 .^ ( -25 / 20 ) * ones(freqbins,N); % 初始化有语音时的条件增益:幅度
Gmin = -30;   % 无语音时的最小增益值:dB
Gainf = 10 ^ ( Gmin / 20 ) * ones(freqbins,N);




tmp_bkg_sp = zeros(freqbins,N);
smoothed_tmp_bkg_sp = zeros(freqbins,N);
tmp_Ef_Prev = zeros(freqbins,N);
min_Ef_Prev = zeros(freqbins,N);
NoiseSpectrumEst_Prev = zeros(freqbins,N);
InstantSmoothedPf = zeros(freqbins,N);
RecursiveEf =  zeros(freqbins,N);
RecursiveEf_Prev =  zeros(freqbins,N);
ratio =  zeros(freqbins,N);
SpeechPresenceProb_Prev = 0.01 * ones(freqbins,N);
RecurAvgPriorSNR =  zeros(freqbins,N);
PriorSNR_Prev = 10 .^ ( -25 / 20 ) * ones(freqbins, N);
RecurAvgPriorSNR_PeakPrev = 10 ^ ( -25 / 20 ) * ones(1, N);
AvgRecurAvgPriorSNRPrev   = 10 ^ ( -25 / 20 ) * ones(1, N);
RecurAvgPriorSNRPrev = 10 .^ ( -25 / 20 ) * ones(freqbins, N);
mu_l = zeros(1, N);
P_Frame = zeros(1, N);
RecurAvgPriorSNR_Peak = zeros(1, N);
SpeechAbsenceProb = zeros(freqbins,N);
NoiseSpectrumEst = zeros(freqbins,N);
pkl_prev = 0.05 * ones(freqbins,1);
pkl_prev_course = [];
alpha_pkl = 0.8;
%%  vad 参数
f = ( 0 : block_size_ / 2 ) * fs / block_size_;
frequency_start = 50;
frequency_end = 4000;
alpha_log = 2;
frequency_start_ind = ceil(frequency_start * block_size_ / fs);
frequency_end_ind = ceil(frequency_end * block_size_ / fs);
beyond_ratio_energy = 1.2;
beyond_ratio_ee = 1;
avg_pkl_seg_th = 0.2;
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
avg_pkl_seg_prev = 0;

%% doa 参数
theta = ( 0 : 10 : 180 ) * pi / 180;
Lt = length(theta);
freq = ( 0 : block_size_ / 2 ) * fs / block_size_;
frequency_start_doa = 50;
frequency_end_doa = 4000;
frequency_start_ind_doa = ceil(frequency_start_doa * block_size_ / fs);
frequency_end_ind_doa = ceil(frequency_end_doa * block_size_ / fs);
centered = mean(array_shape,1);
centered_array = array_shape - repmat(centered,N,1);
sound_path = zeros(N,1);
sound_speed = 340;

fft_in_block_prev = zeros(freqbins,N);
alpha_doa = 0.8;
%% bf参数
bf_fft_in_block = zeros(freqbins,N);
bf_fft_out_block = zeros(freqbins,N);


%% pf参数
Nc = N * ( N - 1 ) / 2;
auto_sp_prev = zeros(freqbins,N);
inter_sp_prev = zeros(freqbins,Nc);
auto_sp = zeros(freqbins,N);
inter_sp = zeros(freqbins,Nc);
alpha_pf = 0.8; % 递归平滑因子
final_output_buffer_ = zeros( element_count, 1 );
final_out_buf_course = [];
f1 = ceil( block_size_ / fs * 200 );
f2 = ceil( block_size_ / fs * 6800 );
fspan = [ f1 : 1 : f2 ].';
H_min_func = [ zeros( f1 - 1, 1);...
    0.7 / ( f2 - f1 ) .* ( fspan - f1 ); ...
    0.7 .* ones( freqbins - f2, 1 ) ]; % see diploma thesis Boigner






%% 变量历程
InstantSmoothedPf_course = [];
RecursiveEf_Prev_course = [];
bkg_noise_energy_course = [];
bkg_noise_energy_entropy_ratio_course = [];
instant_energy_course = [];
instant_energy_entropy_ratio_course = [];
avg_pkl_seg_course = [];
is_speech_frame_course = [];
is_speech_triggered_course = [];
in_fd_block_course = [];
theta_est_course = zeros(frame_num,1);


%% process
out_buf_course = [];
tic
for i = 1 : frame_num
    is_speech_frame = 0;
    
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind, : );
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
        input_block_ = input_block_ .* repmat(wint,1,N);
        fft_in_block = fft(input_block_);
        fft_in_block = fft_in_block(1:freqbins,:);
        bf_fft_in_block = fft_in_block;
        in_fd_block = abs(fft_in_block(1:freqbins,:)) .^ 2;
        in_fd_block_course = [in_fd_block_course in_fd_block(:,1)];
        ang_in_fd_block = angle(fft_in_block(1:freqbins,:));
        if ( i <= NIS )
            tmp_bkg_sp = tmp_bkg_sp + in_fd_block;
            blocks = blocks + 1;
            % 认为前导帧均是噪声帧
            is_speech_frame = 0;
            is_speech_triggered = 0;
            bf_fft_out_block = bf_fft_in_block;
        else
            %% ns
            if ( i == NIS + 1 )
                tmp_bkg_sp = tmp_bkg_sp / blocks;
                smoothed_tmp_bkg_sp = smoothPf( tmp_bkg_sp, omega,...
                    pf_smooth_win , freqbins, N);
                RecursiveEf_Prev = smoothed_tmp_bkg_sp;
                tmp_Ef_Prev = smoothed_tmp_bkg_sp;
                min_Ef_Prev = smoothed_tmp_bkg_sp;
                NoiseSpectrumEst_Prev = tmp_bkg_sp;
            end
            InstantSmoothedPf = smoothPf( in_fd_block, omega,...
                pf_smooth_win , freqbins, N);
            InstantSmoothedPf_course = [InstantSmoothedPf_course InstantSmoothedPf(:,1)];
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
            RecursiveEf_Prev_course = [RecursiveEf_Prev_course RecursiveEf_Prev(:,1)];
            min_Ef_Prev            = min_Ef;
            tmp_Ef_Prev            = tmp_Ef;
            SpeechPresenceProb_Prev = SpeechPresenceProb;
            NoiseSpectrumEst_Prev   = NoiseSpectrumEst;
            RecurAvgPriorSNR = beta .* RecurAvgPriorSNRPrev + ...
                ( 1 - beta ) .* PriorSNR_Prev;
            % local
            RecurAvgPriorSNR_Local = ...
                smoothPf( RecurAvgPriorSNR, omega_local, h_local , freqbins, N);
            P_Local = RecurAvgPriorSNR_Local;
            % global
            RecurAvgPriorSNR_Global = ...
                smoothPf( RecurAvgPriorSNR, omega_global, h_global , freqbins, N);
            P_Global = RecurAvgPriorSNR_Global;
            for j = 1 : N
                for k = 1 : freqbins
                    if( P_Local( k, j ) <= 10 ^ ( sigma_min / 20 ) )
                        P_Local( k, j ) = 0;
                    elseif( P_Local( k, j ) >= 10 ^ ( sigma_max / 20 ) )
                        P_Local( k, j ) = 1;
                    else
                        P_Local( k, j ) = log10( P_Local( k, j ) / 10 ^ ( sigma_min / 20 ) ) ...
                            / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
                    end
                    
                    if( P_Global( k, j ) <= 10 ^ ( sigma_min / 20 ) )
                        P_Global( k, j ) = 0;
                    elseif( P_Global( k, j ) >= 10 ^ ( sigma_max / 20 ) )
                        P_Global( k, j ) = 1;
                    else
                        P_Global( k, j ) = log10( P_Global( k, j ) / 10 ^ ( sigma_min / 20 ) ) ...
                            / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
                    end
                end
            end
            AvgRecurAvgPriorSNR = mean( RecurAvgPriorSNR( 1 : ( frame_shift / 2 + 1), : ), 1 );
            th1 = 10 ^ ( sigma_min / 20 ) * 10 .^ ( RecurAvgPriorSNR_PeakPrev / 20 );
            th2 = 10 ^ ( sigma_max / 20 ) * 10 .^ ( RecurAvgPriorSNR_PeakPrev / 20 );
            % 计算mu_l
            for j = 1 : N
                if( AvgRecurAvgPriorSNR(j) <= th1(j) )
                    mu_l(j) = 0;
                elseif( AvgRecurAvgPriorSNR(j) >= th2(j) )
                    mu_l(j) = 1;
                else
                    mu_l(j) = log10( AvgRecurAvgPriorSNR(j) / th1(j) ) ...
                        / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
                end
            end
            % 计算P_Frame
            for j = 1 : N
                if( AvgRecurAvgPriorSNR(j) <= 10 ^ ( sigma_min / 20 ) )
                    P_Frame(j) = 0;
                elseif( AvgRecurAvgPriorSNR(j) <= AvgRecurAvgPriorSNRPrev(j) )
                    P_Frame(j) = mu_l(j);
                else
                    RecurAvgPriorSNR_Peak(j) = min( max( AvgRecurAvgPriorSNR(j), 10 ^ ( sigma_pmin / 20 ) ),...
                        10 ^ ( sigma_pmax / 20 ) );
                    RecurAvgPriorSNR_PeakPrev(j) = RecurAvgPriorSNR_Peak(j);
                    P_Frame(j) = 1;
                end
            end
            % 计算无语音概率
            SpeechAbsenceProb = 1 - ...
                P_Local .* P_Global .* repmat(P_Frame,freqbins,1);
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
            full_com_sp = [half_com_sp;conj(half_com_sp(end-1:-1:2,:))];
            output_block_ = real(ifft(full_com_sp));
            
            %% vad
            % 计算两个通道噪声谱和瞬时信号谱的和,以此进行语音检测
            chs_avg_sf = sum(RecursiveEf,2);
%             chs_avg_sf = sum(InstantSmoothedPf,2);
            chs_avg_nf = sum(NoiseSpectrumEst,2);
%             chs_avg_nf = sum(smoothed_tmp_bkg_sp,2);
            % 计算背景噪声以及含噪信号的能量和能熵比判决统计量
            % 背景噪声
            bkg_noise_energy = mean(chs_avg_nf...
                (frequency_start_ind:frequency_end_ind));
            bkg_noise_energy_entropy_ratio = ...
                ComputeLogEnergyEntropyRatio(chs_avg_nf,...
                frequency_start_ind,frequency_end_ind,alpha_log);
            bkg_noise_energy_course = [bkg_noise_energy_course bkg_noise_energy];
            bkg_noise_energy_entropy_ratio_course = [bkg_noise_energy_entropy_ratio_course...
                bkg_noise_energy_entropy_ratio];
            % 含噪信号
            instant_energy = mean(chs_avg_sf...
                (frequency_start_ind:frequency_end_ind));
            instant_energy_entropy_ratio = ComputeLogEnergyEntropyRatio(chs_avg_sf,...
                frequency_start_ind,frequency_end_ind,alpha_log);
            instant_energy_course = [instant_energy_course instant_energy];
            instant_energy_entropy_ratio_course = [instant_energy_entropy_ratio_course ...
                instant_energy_entropy_ratio];
             % 判断该帧是不是语音帧
             % 根据能量
%             if (instant_energy > beyond_ratio_energy * bkg_noise_energy)
%                 is_speech_frame_block = 1;
%             else
%                 is_speech_frame_block = 0;
%             end
            % 根据能熵比
%             if (instant_energy_entropy_ratio > beyond_ratio_ee * bkg_noise_energy_entropy_ratio)
%                 is_speech_frame_block = 1;
%             else
%                 is_speech_frame_block = 0;
%             end
            % 根据语音存在概率
            avg_pkl = mean(pkl,2);
            pkl_smooth = smoothPf( avg_pkl, omega_local, h_local , freqbins, 1);
            
            pkl_smooth = alpha_pkl * pkl_prev + ( 1 - alpha_pkl ) * pkl_smooth;
            pkl_prev = pkl_smooth;
            pkl_prev_course = [pkl_prev_course  pkl_smooth];
            sort_pkl = sort(pkl_smooth,'descend');
            avg_pkl_seg = alpha_pkl * avg_pkl_seg_prev + ... 
                ( 1 - alpha_pkl ) * mean(mean(sort_pkl(1:10,:)));
            avg_pkl_seg_prev = avg_pkl_seg;
            
            if (avg_pkl_seg > avg_pkl_seg_th)
                is_speech_frame_block = 1;
            else
                is_speech_frame_block = 0;
            end
            avg_pkl_seg_course = [avg_pkl_seg_course avg_pkl_seg];
            is_speech_frame = is_speech_frame || is_speech_frame_block;
            
            %% doa 
            fft_in_block = alpha_doa * fft_in_block_prev + ...
                ( 1 - alpha_doa ) * fft_in_block .* Gainf ;
            fft_in_block_prev = fft_in_block;
            if (is_speech_frame)
                Pf_out = zeros(Lt,1);
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
                    
                    for m = frequency_start_ind_doa : frequency_end_ind_doa;
                        fm = freq(m);
                        sv = exp( -1j * 2 * pi * fm * tao );
                        xfk = fft_in_block(m,:).'; % 单频点谱值向量
                        %% cbf
                        Pf_out(j) = Pf_out(j) + abs(sv' * ( xfk * xfk' ) * sv);
                        %% mvdr
%                         Rxf_mvdr = xfk * xfk' + 0.01 * eye( N );
%                         Pf_out(j) = Pf_out(j) + 1 / abs( sv' / Rxf_mvdr * sv );
                        %% stmv
%                         Tf = diag(sv);
%                         Rxf_stmv = xfk * xfk';
%                         RSCTM = RSCTM + Tf' * Rxf_stmv * Tf;
                    end  % end of frequency scan 
%                     RSCTM = RSCTM + 0.01 * eye( N );
%                     wo = ( RSCTM \ I ) / ( I' / RSCTM * I );
%                     Pf_out(j) = 1 / abs( I' / RSCTM * I );
                    
                end  % end of Lt
                [max_pf,max_index] = max(Pf_out);
                theta_est = theta(max_index) * 180 / pi;
                theta_est_course(i) = theta_est;
              
                
                %% bf
                %% 常规波束形成
                for k = 1 : N
                    sound_path(k) = centered_array(k,1) * cosd(theta_est) + ...
                        centered_array(k,2) * sind(theta_est);
                end
                tao = sound_path / sound_speed;
                for j = 1 : freqbins
                    fm = freq(j);
                    sv = exp( -1j * 2 * pi * fm * ( tao - tao(1) ) );
                    xfk = bf_fft_in_block(j,:).'; % 单频点谱值向量
                    bf_fft_out_block(j,:) = conj(sv) .* xfk;
                end
            else
                bf_fft_out_block = bf_fft_in_block;
            end  % end of if is_speech_frame
        end  % end of if i <= NIS
        
        
        bf_fft_out_block = bf_fft_out_block .* Gainf;
        %% 后置滤波
        %% 自功率谱估计
        auto_sp = bf_fft_out_block .* conj(bf_fft_out_block);
        auto_sp = alpha_pf * auto_sp_prev + ( 1 - alpha_pf ) * auto_sp;
        auto_sp_prev = auto_sp;
        
        %% 互功率谱估计
        k = 1;
        for k1 = 1 : N - 1
            for k2 = k1 + 1 : N
                inter_sp(:,k) = real( bf_fft_out_block(:,k1) .* ... 
                    conj(bf_fft_out_block(:,k2)));
                k = k + 1;
            end
        end
        inter_sp = alpha_pf * inter_sp_prev + ( 1 - alpha_pf ) * inter_sp;
        inter_sp_prev = inter_sp;
        
        %% 维纳滤波
        Ap = mean(auto_sp,2);
        Cp = mean(inter_sp,2);
        Hz = Cp ./ Ap;
        Hz = max(Hz,0.05);
%         Hz = max(Hz,H_min_func);
        Hz = min(Hz,1);
        
        %% 直接平均多通道频谱(包括相位谱)的方式
%         final_fft_out_block = mean(bf_fft_out_block,2);
%         final_fft_out_block = Hz .* final_fft_out_block;
        %% 平均多通道能量谱，利用某一通道的相位谱
        square_sp = abs(bf_fft_out_block) .^ 2;
        avg_square_sp = sum(square_sp,2);
        angle_sp = angle(bf_fft_out_block);
        final_fft_out_block = Hz .* sqrt(avg_square_sp) ...
            .* exp( 1j * angle_sp(:,1) );
        
        %% 构造复数全谱
        full_com_sp = [final_fft_out_block;conj(final_fft_out_block(end-1:-1:2,:))];
        final_output_block_ = real(ifft(full_com_sp,[],1));
         %% 加窗
        final_output_block_ = final_output_block_ .* wint;
        
        
        
        output_block_ = output_block_ .* repmat(wint,1,N);
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
        
        final_output_buffer_(result_ind_range, :) =  ...
            final_output_buffer_(a_ind_range, :) + ...
            final_output_block_(b_ind_range,:);
        
        first_frame_in_block = first_frame_in_block + shift_amount_;
    end  % end of while
    
    
    
    % 判断该帧是否在有效语音段内
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
    
    is_speech_frame_course = [is_speech_frame_course is_speech_frame];
    is_speech_triggered_course = [is_speech_triggered_course is_speech_triggered];
    
    
    % Copy output buffer to output
    % CopyFrames
    out_buf = output_buffer_(1:chunk_length_,:);
    final_out_buf = final_output_buffer_(1:chunk_length_,:);
    % MoveFrames ----> output_buffer
    dst_ind_range = 1 : initial_delay_;
    src_ind_range = ( chunk_length_ + 1 ) : ...
        ( chunk_length_ + initial_delay_ );
    output_buffer_(dst_ind_range,:) = output_buffer_(src_ind_range,:);
    final_output_buffer_(dst_ind_range,:) = final_output_buffer_(src_ind_range,:);
    % ZeroOut ----> output_buffer
    src_ind_range = ( initial_delay_ + 1 ) : ...
        ( chunk_length_ + initial_delay_ );
    output_buffer_(src_ind_range,:) = 0;
    final_output_buffer_(src_ind_range,:) = 0;
    frame_offset_ = first_frame_in_block - chunk_length_;
    out_buf_course = [out_buf_course;out_buf];
    final_out_buf_course = [final_out_buf_course;final_out_buf];
    
end
toc

%% 绘制降噪波形
LL = length(out_buf_course(225:end,1));
figure
plot((1:L)/fs,audio(:,1))
hold on
plot((1:LL)/fs,out_buf_course(225:end,1),'r')
grid on
axis tight
legend('模拟输入序列','未经算法处理的输出序列')

%% 绘制最终输出的单通道波形
Lout = length(final_out_buf_course);
figure
plot((1:Lout)/fs,final_out_buf_course)
grid on
axis tight
title('阵列信号处理输出')

audiowrite('array_processed_out.wav',final_out_buf_course,fs);

%% 绘制vad相关统计量变化曲线
figure
subplot 211
plot((1:frame_num)*0.01, is_speech_frame_course)
axis tight
ylim([-1 2])
subplot 212
plot((1:frame_num)*0.01, is_speech_triggered_course)
axis tight
ylim([-1 2])

figure
plot(avg_pkl_seg_course)
axis tight
title('语音存在概率')

figure
[M1,M2] = size(pkl_prev_course);
imagesc(1:M2,f,pkl_prev_course);
axis xy
xlabel('帧序号/n')
ylabel('f/Hz')

figure
imagesc(1:M2,f,10*log10(in_fd_block_course));
% imagesc(1:M2,f,10*log10(InstantSmoothedPf_course));
% imagesc(1:M2,f,10*log10(RecursiveEf_Prev_course));
axis xy
xlabel('帧序号/n')
ylabel('f/Hz')



figure
subplot 211
plot(instant_energy_course)
hold on
plot(bkg_noise_energy_course,'r')
plot(bkg_noise_energy_course * beyond_ratio_energy ,'g')
axis tight
title('瞬时能量和')
subplot 212
plot(bkg_noise_energy_course)
axis tight
title('参考背景能量和')

figure
subplot 211
plot(instant_energy_entropy_ratio_course)
hold on
plot(bkg_noise_energy_entropy_ratio_course,'r')
plot(bkg_noise_energy_entropy_ratio_course * beyond_ratio_ee ,'g')
axis tight
title('瞬时能熵比')
subplot 212
plot(bkg_noise_energy_entropy_ratio_course)
axis tight
title('参考背景能熵比')

audiowrite('ns.wav',out_buf_course,fs);

figure
plot(1:frame_num,theta_est_course);
axis tight
xlabel('帧序号/n')
ylabel('方位估计值/°')





