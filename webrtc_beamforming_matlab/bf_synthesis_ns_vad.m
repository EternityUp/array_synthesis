% ���������롢�������⡢��Դ���򡢲����γɺͺ����˲���ϵ�һ��
clear all
close all
clc
%% ��ȡ��������
[audio,fs] = audioread('stereo.wav');
% audio = audio(:,1);
[M,N] = size(audio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_time = 10; % ÿ֡ʱ��
frame_len = fs * frame_time / 1000; % ÿ֡����
frame_shift = 128; % ����֮֡���֡��
lapped_len = frame_len - frame_shift; % ����֡���ص�����
% nfft = 2 ^ nextpow2( frame_len ); % ����Ҷ�任����
nfft = 2 * frame_len - 1;
% freq_num = 1 + nfft / 2; % ����Ƶ�ʵ���
freq_num = frame_len;
wint = hanning( frame_len ); % ������
f = ( 0 : nfft / 2 ) * fs / nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_num = 0; % ��ʼ��֡��
audio_ret = 0; % ��ʼ����ʼ��ָ��λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ȹ�����س���
NIS = 5; % ǰ������֡��Ŀ�����ڱ��������׹���
bkgPf = zeros( freq_num, N ); % ����������
omega = 3; % Ƶ��ƽ�����ڴ�С,ȡ����
pf_smooth_win = hanning( omega );
% ������������һ��
pf_smooth_win = pf_smooth_win ./ sqrt( pf_smooth_win' * pf_smooth_win );
alpha_s = 0.8; % ������ƽ������
NL = 125; % ���ƺ����ź�frequency-wise������Сֵ��֡����С
delta = 5; % �������Ƶ���Ȼ����ֵ
alpha_p = 0.8; % �����������������ʵ�ƽ������
alpha_d = 0.95 * ones(freq_num, N); % �����׹��Ƶ�ƽ������
beta = 0.7;        % ���Ƶݹ�ƽ����������ȵ�ƽ������
SpeechPresenceProb_Prev = 0.01 * ones(freq_num, N); % ��ʼ��ǰһ֡����������������
PriorSNR_Prev = 10 .^ ( -25 / 20 ) * ones(freq_num, N);         % ��ʼ�����������:����
RecurAvgPriorSNRPrev = 10 .^ ( -25 / 20 ) * ones(freq_num, N);  % ��ʼ���ݹ�ƽ�������������:����
RecurAvgPriorSNR_PeakPrev = 10 ^ ( -25 / 20 ) * ones(1, N);
AvgRecurAvgPriorSNRPrev   = 10 ^ ( -25 / 20 ) * ones(1, N); % ĳƵ���ڵĵݹ�ƽ������������Ⱦ�ֵ:����
omega_local = 3;   % �ֲ�ƽ�����ڴ�С
h_local = hanning(omega_local);  % �ֲ�ƽ��������
h_local = h_local ./ sqrt( h_local' * h_local ); % ������������һ��
omega_global = 31; % ȫ��ƽ�����ڴ�С
h_global = hanning(omega_global); % ȫ��ƽ��������
h_global = h_global ./ sqrt( h_global' * h_global ); % ������������һ��
sigma_min = -10;   % ˥������ͬʱ������������������������Ⱦ��鳣��:dB
sigma_max = -5;    % ˥������ͬʱ������������������������Ⱦ��鳣��:dB
sigma_pmin = 0;    % ���ڿ�����������������ת�Ƶ���������Ⱦ��鳣��:dB
sigma_pmax = 10;   % ���ڿ�����������������ת�Ƶ���������Ⱦ��鳣��:dB
q_max = 0.95;      % ������������ʵ����ֵ
alpha = 0.99 * ones(freq_num,N); % ����ƽ�����������̶�����ɽ׶����������ʧ���ƽ������
PostSNR_Prev  = 10 .^ ( -25 / 20 ) * ones(freq_num,N); % ��ʼ�����������:����
GainH1_Prev   = 10 .^ ( -25 / 20 ) * ones(freq_num,N); % ��ʼ��������ʱ����������:����
pkl_prev = zeros(freq_num,N);
Gmin = -30;   % ������ʱ����С����ֵ:dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAD��ز���
% ѡ��ĳ������Ƶ�Σ������������ʵ�ƽ��ֵ
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
        xin = audio( iniind : endind, : ); % ��ȡÿ֡����
    else % ���һ�����ݲ���һ֡
        endind = M;
        zero_num = frame_len - endind + iniind - 1;
        xin = [audio( iniind : endind, : ); ...
            zeros(zero_num, N)]; % ��ȡÿ֡����
    end
    frame_num = frame_num + 1;
    xin_win = xin .* repmat( wint, 1, N );
    %% ʱƵ��ת��
    fft_xin_win = fft( xin_win, nfft, 1);
    pf_xin = fft_xin_win( 1 : freq_num, : );
    %% �ж��ǲ���ǰ��֡����������Ӧ����
    if ( frame_num <= NIS )
        bkgPf = bkgPf + abs(pf_xin) .^ 2 / NIS;
    end
    %% ��������Ƶ��ƽ��
    if ( frame_num == NIS )
        smoothedPf = smoothPf( bkgPf, omega, pf_smooth_win , freq_num, N);
        RecursiveEf_Prev = smoothedPf;
        tmp_Ef_Prev = smoothedPf;
        min_Ef_Prev = smoothedPf;
        NoiseSpectrumEst_Prev = bkgPf; % ��ʼ��ǰһ֡�����������׹���
    end
    
    if ( frame_num > NIS )
        % �Ӵ�ƽ����ǰ֡�źŹ�����
        Ef = abs(pf_xin) .^ 2;
        smoothedEf = smoothPf( Ef, omega, pf_smooth_win , freq_num, N);
        % ����ݹ�ƽ�����������ס���Сֵ��˲ʱֵ
        RecursiveEf = alpha_s * RecursiveEf_Prev + ...
            ( 1 - alpha_s ) * smoothedEf;
        min_Ef = min( min_Ef_Prev, RecursiveEf );
        tmp_Ef = min( tmp_Ef_Prev, RecursiveEf );
        if ( mod( frame_num, NL ) == 0 )
            min_Ef = min( tmp_Ef_Prev, RecursiveEf );
            tmp_Ef = RecursiveEf;
        end
        % ����RecursiveEf��min_Ef�ı�ֵ
        ratio = RecursiveEf ./ min_Ef;
        % ��������������
        Ikl = ratio > delta;
        SpeechPresenceProb = alpha_p .* SpeechPresenceProb_Prev + ...
            ( 1 - alpha_p ) * Ikl;
        % �������������׹��Ƶ�ʱ��ƽ������
        alpha_d = alpha_d + ( 1 - alpha_d ) .* SpeechPresenceProb;
        NoiseSpectrumEst = alpha_d .* NoiseSpectrumEst_Prev + ...
            ( 1 - alpha_d ) .* smoothedEf;
        % --------����ǰһ֡�����������׹�����ر���------------------------%
        RecursiveEf_Prev        = RecursiveEf;
        min_Ef_Prev            = min_Ef;
        tmp_Ef_Prev            = tmp_Ef;
        SpeechPresenceProb_Prev = SpeechPresenceProb;
        NoiseSpectrumEst_Prev   = NoiseSpectrumEst;
        % ��ǰ֡�������������ʹ���
        % ���㵱ǰ֡��������ȵĵݹ�ƽ��
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
            % ����mu_l
            if( AvgRecurAvgPriorSNR(j) <= th1(j) )
                mu_l(j) = 0;
            elseif( AvgRecurAvgPriorSNR(j) >= th2(j) )
                mu_l(j) = 1;
            else
                mu_l(j) = log10( AvgRecurAvgPriorSNR(j) / th1(j) ) ...
                    / log10( 10 ^ ( sigma_max / 20 ) / 10 ^ ( sigma_min / 20 ) );
            end
            % ����P_Frame
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
            % ��������������
            SpeechAbsenceProb(:,j) = 1 - ... 
                P_Local(:,j) .* P_Global(:,j) .* P_Frame(j);
        end
        % �����������������ֵΪ0.95
        SpeechAbsenceProb( SpeechAbsenceProb >= q_max ) = q_max;

        % ------------������ر���-------------------%
        AvgRecurAvgPriorSNRPrev = AvgRecurAvgPriorSNR;
        RecurAvgPriorSNRPrev = RecurAvgPriorSNR;
        % ������������
        PostSNRCurr = Ef ./ NoiseSpectrumEst; % ����
        
        % ����ǰһ֡���������������桢��������Ⱥ͵�ǰ֡�ĺ��������,
        % ���Ƶ�ǰ֡�����������
        PriorSNRCurr = alpha .* ( GainH1_Prev .^ 2 ) .* PostSNR_Prev + ...
            ( 1 - alpha ) .* max( PostSNRCurr - 1, 0 );
        vkPost_Prior = PriorSNRCurr .* PostSNRCurr ./ ( 1 + PriorSNRCurr );
        
        vkPost_Prior = max(vkPost_Prior,1e-11);
        % ���Ƶ�ǰ֡����������������
        GainH1_Curr = PriorSNRCurr ./ ( 1 + PriorSNRCurr ) .* ...
            exp( expint( vkPost_Prior ) / 2 );
        % --------------��������ͺ�������ȵ�----------------------- %
        PriorSNR_Prev = PriorSNRCurr;
        PostSNR_Prev  = PostSNRCurr;
%         GainH1_Prev   = GainH1_Curr;
        GainH1_Prev   = max(GainH1_Curr,10 .^ ( Gmin / 20 ));
     
        
        %% ����������ֵ���м����������ƺ��������
        % ����ԭʼ�����ź��е�����������
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