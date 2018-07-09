% 将语音降噪、声音活动检测、声源测向、波束形成和后置滤波结合到一起
clear all
close all
clc
%% 读取语音数据
[audio,fs] = audioread('stereo.wav');
% audio = audio(:,1);
[M,N] = size(audio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_time = 10; % 每帧时长
frame_len = fs * frame_time / 1000; % 每帧长度
frame_shift = 128; % 连续帧之间的帧移
lapped_len = frame_len - frame_shift; % 连续帧间重叠点数
% nfft = 2 ^ nextpow2( frame_len ); % 傅里叶变换点数
nfft = 2 * frame_len - 1;
% freq_num = 1 + nfft / 2; % 半轴频率点数
freq_num = frame_len;
wint = hanning( frame_len ); % 窗函数
f = ( 0 : nfft / 2 ) * fs / nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_num = 0; % 初始化帧数
audio_ret = 0; % 初始化起始读指针位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 信噪比估计相关常数
NIS = 5; % 前导噪声帧数目，用于背景噪声谱估计
bkgPf = zeros( freq_num, N ); % 背景噪声谱
omega = 3; % 频域平滑窗口大小,取奇数
pf_smooth_win = hanning( omega );
% 窗函数能量归一化
pf_smooth_win = pf_smooth_win ./ sqrt( pf_smooth_win' * pf_smooth_win );
alpha_s = 0.8; % 功率谱平滑因子
NL = 125; % 估计含噪信号frequency-wise能量最小值的帧窗大小
delta = 5; % 噪声估计的似然比阈值
alpha_p = 0.8; % 估计有语音条件概率的平滑因子
alpha_d = 0.95 * ones(freq_num, N); % 噪声谱估计的平滑因子
beta = 0.7;        % 估计递归平均先验信噪比的平滑因子
SpeechPresenceProb_Prev = 0.01 * ones(freq_num, N); % 初始化前一帧的有语音条件概率
PriorSNR_Prev = 10 .^ ( -25 / 20 ) * ones(freq_num, N);         % 初始化先验信噪比:幅度
RecurAvgPriorSNRPrev = 10 .^ ( -25 / 20 ) * ones(freq_num, N);  % 初始化递归平均的先验信噪比:幅度
RecurAvgPriorSNR_PeakPrev = 10 ^ ( -25 / 20 ) * ones(1, N);
AvgRecurAvgPriorSNRPrev   = 10 ^ ( -25 / 20 ) * ones(1, N); % 某频段内的递归平均的先验信噪比均值:幅度
omega_local = 3;   % 局部平滑窗口大小
h_local = hanning(omega_local);  % 局部平滑汉宁窗
h_local = h_local ./ sqrt( h_local' * h_local ); % 窗函数能量归一化
omega_global = 31; % 全局平滑窗口大小
h_global = hanning(omega_global); % 全局平滑汉宁窗
h_global = h_global ./ sqrt( h_global' * h_global ); % 窗函数能量归一化
sigma_min = -10;   % 衰减噪声同时保持弱语音分量的先验信噪比经验常数:dB
sigma_max = -5;    % 衰减噪声同时保持弱语音分量的先验信噪比经验常数:dB
sigma_pmin = 0;    % 用于控制语音段向噪声段转移的先验信噪比经验常数:dB
sigma_pmax = 10;   % 用于控制语音段向噪声段转移的先验信噪比经验常数:dB
q_max = 0.95;      % 无语音先验概率的最大值
alpha = 0.99 * ones(freq_num,N); % 用于平衡噪声消除程度与过渡阶段引入的语音失真的平滑因子
PostSNR_Prev  = 10 .^ ( -25 / 20 ) * ones(freq_num,N); % 初始化后验信噪比:幅度
GainH1_Prev   = 10 .^ ( -25 / 20 ) * ones(freq_num,N); % 初始化有语音时的条件增益:幅度
pkl_prev = zeros(freq_num,N);
Gmin = -30;   % 无语音时的最小增益值:dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAD相关参数
% 选定某个语音频段，计算语音概率的平均值
freq_range = 3 : 17;
speech_presence_threshold = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_l = zeros(N,1);
P_Frame = zeros(N,1);
RecurAvgPriorSNR_Peak = zeros(N,1);
SpeechAbsenceProb = zeros(freq_num,N);
tic
while ( audio_ret < M )
    iniind = audio_ret + 1;
    if ( audio_ret + frame_len <= M )
        endind = iniind + frame_len - 1;
        xin = audio( iniind : endind, : ); % 提取每帧数据
    else % 最后一段数据不足一帧
        endind = M;
        zero_num = frame_len - endind + iniind - 1;
        xin = [audio( iniind : endind, : ); ...
            zeros(zero_num, N)]; % 提取每帧数据
    end
    frame_num = frame_num + 1;
    xin_win = xin .* repmat( wint, 1, N );
    %% 时频域转换
    fft_xin_win = fft( xin_win, nfft, 1);
    pf_xin = fft_xin_win( 1 : freq_num, : );
    %% 判断是不是前导帧，并进行相应处理
    if ( frame_num <= NIS )
        bkgPf = bkgPf + abs(pf_xin) .^ 2 / NIS;
    end
    %% 背景噪声频谱平滑
    if ( frame_num == NIS )
        smoothedPf = smoothPf( bkgPf, omega, pf_smooth_win , freq_num, N);
        RecursiveEf_Prev = smoothedPf;
        tmp_Ef_Prev = smoothedPf;
        min_Ef_Prev = smoothedPf;
        NoiseSpectrumEst_Prev = bkgPf; % 初始化前一帧的噪声能量谱估计
    end
    
    if ( frame_num > NIS )
        % 加窗平滑当前帧信号功率谱
        Ef = abs(pf_xin) .^ 2;
        smoothedEf = smoothPf( Ef, omega, pf_smooth_win , freq_num, N);
        % 计算递归平滑噪声能量谱、最小值、瞬时值
        RecursiveEf = alpha_s * RecursiveEf_Prev + ...
            ( 1 - alpha_s ) * smoothedEf;
        min_Ef = min( min_Ef_Prev, RecursiveEf );
        tmp_Ef = min( tmp_Ef_Prev, RecursiveEf );
        if ( mod( frame_num, NL ) == 0 )
            min_Ef = min( tmp_Ef_Prev, RecursiveEf );
            tmp_Ef = RecursiveEf;
        end
        % 计算RecursiveEf与min_Ef的比值
        ratio = RecursiveEf ./ min_Ef;
        % 估计有语音概率
        Ikl = ratio > delta;
        SpeechPresenceProb = alpha_p .* SpeechPresenceProb_Prev + ...
            ( 1 - alpha_p ) * Ikl;
        % 更新噪声能量谱估计的时变平滑参数
        alpha_d = alpha_d + ( 1 - alpha_d ) .* SpeechPresenceProb;
        NoiseSpectrumEst = alpha_d .* NoiseSpectrumEst_Prev + ...
            ( 1 - alpha_d ) .* smoothedEf;
        % --------更新前一帧的噪声能量谱估计相关变量------------------------%
        RecursiveEf_Prev        = RecursiveEf;
        min_Ef_Prev            = min_Ef;
        tmp_Ef_Prev            = tmp_Ef;
        SpeechPresenceProb_Prev = SpeechPresenceProb;
        NoiseSpectrumEst_Prev   = NoiseSpectrumEst;
        % 当前帧先验无语音概率估计
        % 计算当前帧先验信噪比的递归平均
        RecurAvgPriorSNR = beta .* RecurAvgPriorSNRPrev + ...
            ( 1 - beta ) .* PriorSNR_Prev;
        
        % local 
        RecurAvgPriorSNR_Local = ...
            smoothPf( RecurAvgPriorSNR, omega_local, h_local , freq_num, N);
        P_Local = RecurAvgPriorSNR_Local;
        % global
        RecurAvgPriorSNR_Global = ...
            smoothPf( RecurAvgPriorSNR, omega_global, h_global , freq_num, N);
        P_Global = RecurAvgPriorSNR_Global;
        
        for j = 1 : N
            for k = 1 : freq_num
                if( P_Local( k, j ) <= 10 ^ ( sigma_min / 20 ) )
                    P_Local( k, j ) = 0;
                elseif( P_Local( k, j ) >= 10 ^ ( sigma_max / 20 ) )
                    P_Local( k, j ) = 1;
                else
                    P_Local( k, j ) = log10( P_Local( k, j ) / 10 ^ ( sigma_min / 20 ) ) ...
                        / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
                end
            end
            
            for k = 1 : freq_num
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
        
        AvgRecurAvgPriorSNR = mean( RecurAvgPriorSNR( 1 : ( frame_shift / 2 ), : ), 1 );
        th1 = 10 ^ ( sigma_min / 20 ) * 10 .^ ( RecurAvgPriorSNR_PeakPrev / 20 );
        th2 = 10 ^ ( sigma_max / 20 ) * 10 .^ ( RecurAvgPriorSNR_PeakPrev / 20 );
        for j = 1 : N
            % 计算mu_l
            if( AvgRecurAvgPriorSNR(j) <= th1(j) )
                mu_l(j) = 0;
            elseif( AvgRecurAvgPriorSNR(j) >= th2(j) )
                mu_l(j) = 1;
            else
                mu_l(j) = log10( AvgRecurAvgPriorSNR(j) / th1(j) ) ...
                    / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
            end
            % 计算P_Frame
            if( AvgRecurAvgPriorSNR(j) < 10 ^ ( sigma_min / 20 ) )
                P_Frame(j) = 0;
            elseif( AvgRecurAvgPriorSNR(j) > AvgRecurAvgPriorSNRPrev(j) )
                P_Frame(j) = mu_l(j);
            else
                RecurAvgPriorSNR_Peak(j) = min( max( AvgRecurAvgPriorSNR(j),10 ^ ( sigma_pmin / 20 ) ),...
                    10 ^ ( sigma_pmax / 20 ) );
                RecurAvgPriorSNR_PeakPrev(j) = RecurAvgPriorSNR_Peak(j);
                P_Frame(j) = 1;
            end
            % 计算无语音概率
            SpeechAbsenceProb(:,j) = 1 - ... 
                P_Local(:,j) .* P_Global(:,j) .* P_Frame(j);
        end
        % 限制无语音概率最大值为0.95
        SpeechAbsenceProb( SpeechAbsenceProb >= q_max ) = q_max;

        % ------------更新相关变量-------------------%
        AvgRecurAvgPriorSNRPrev = AvgRecurAvgPriorSNR;
        RecurAvgPriorSNRPrev = RecurAvgPriorSNR;
        % 计算后验信噪比
        PostSNRCurr = Ef ./ NoiseSpectrumEst; % 幅度
        
        % 根据前一帧的有语音条件增益、后验信噪比和当前帧的后验信噪比,
        % 估计当前帧的先验信噪比
        PriorSNRCurr = alpha .* ( GainH1_Prev .^ 2 ) .* PostSNR_Prev + ...
            ( 1 - alpha ) .* max( PostSNRCurr - 1, 0 );
        vkPost_Prior = PriorSNRCurr .* PostSNRCurr ./ ( 1 + PriorSNRCurr );
        
        vkPost_Prior = max(vkPost_Prior,1e-11);
        % 估计当前帧的有语音条件增益
        GainH1_Curr = PriorSNRCurr ./ ( 1 + PriorSNRCurr ) .* ...
            exp( expint( vkPost_Prior ) / 2 );
        % --------------更新先验和后验信噪比等----------------------- %
        PriorSNR_Prev = PriorSNRCurr;
        PostSNR_Prev  = PostSNRCurr;
%         GainH1_Prev   = GainH1_Curr;
        GainH1_Prev   = max(GainH1_Curr,10 .^ ( Gmin / 20 ));
     
        
        %% 根据谱增益值序列计算噪声抑制后的语音谱
        % 计算原始含噪信号中的有语音概率
        pkl = 1 ./ ( 1 + SpeechAbsenceProb ./ ( 1 - SpeechAbsenceProb ) ...
            .* ( 1 + PriorSNRCurr ) .* exp( -vkPost_Prior ) );
        pkl = 0.8 * pkl + 0.2 * pkl_prev;
        pkl_prev = pkl;
        pkl_course(frame_num,:,:) = pkl;
        Ef_course(frame_num,:,:) = Ef;
        post_snr_course(frame_num,:,:) = PostSNRCurr;
        
        %% simple VAD
        vad_decision(frame_num,:) = mean(pkl(freq_range,:)) > ... 
            speech_presence_threshold;
    end 
    
    audio_ret = audio_ret + frame_shift;
    
    
end
toc
frame_t = ( 0 : ( frame_num - 1 ) ) * frame_shift / fs;
t = ( 1 : M ) / fs;

figure
subplot 211
hold on
plot(frame_t,vad_decision(:,1),'b*')
plot(frame_t,vad_decision(:,2),'ro')
axis([min(frame_t) max(frame_t) -1 2])
legend('ch1','ch2')
xlabel('time/s')
ylabel('decision(0-no;1-yes)')
title('vad decision')
subplot 212
hold on
plot(t,audio(:,1),'b.')
plot(t,audio(:,2),'r--')
legend('ch1','ch2')
xlabel('time/s')
ylabel('amplitude')
xlim([min(t) max(t)])
title('wave')


figure
subplot 211
imagesc(frame_t,f,pkl_course(:,:,1)')
axis xy
xlabel('time/s')
ylabel('frequency/Hz')
title('Channel 1: speech presence probability')
subplot 212
imagesc(frame_t,f,pkl_course(:,:,2)')
axis xy
xlabel('time/s')
ylabel('frequency/Hz')
title('Channel 2: speech presence probability')


Ef_course = 10*log10(abs(Ef_course));
figure
subplot 211
imagesc(frame_t,f,Ef_course(:,:,1)')
axis xy
xlabel('time/s')
ylabel('frequency/Hz')
title('Channel 1: speech presence probability')
subplot 212
imagesc(frame_t,f,Ef_course(:,:,2)')
axis xy
xlabel('time/s')
ylabel('frequency/Hz')
title('Channel 2: speech presence probability')


% figure
% subplot 211
% imagesc(frame_t,f,post_snr_course(:,:,1)')
% axis xy
% xlabel('time/s')
% ylabel('frequency/Hz')
% axis([min(frame_t) max(frame_t) min(f) max(f)])
% 
% subplot 212
% imagesc(frame_t,f,post_snr_course(:,:,2)')
% axis xy
% xlabel('time/s')
% ylabel('frequency/Hz')
% axis([min(frame_t) max(frame_t) min(f) max(f)])
% title('Channel 2: snr')

save vad_decision.mat vad_decision
save post_snr_course.mat post_snr_course