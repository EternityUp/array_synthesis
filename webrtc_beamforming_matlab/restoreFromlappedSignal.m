% �ص��������źŻָ�
clear all
close all
clc
%% ��ȡ��������
[audio,fs] = audioread('stereo.wav');
[M,N] = size(audio);
frame_time = 10; % ÿ֡ʱ��
frame_len = fs * frame_time / 1000; % ÿ֡����
frame_shift = 128; % ����֮֡���֡��
lapped_len = frame_len - frame_shift; % ����֡���ص�����
frame_num = floor( ( M - lapped_len ) / frame_shift ); % ֡��
mono_audio = audio(:,1);
nfft = 2 ^ nextpow2( frame_len );
freq_num = 1 + nfft / 2;
wint = hanning(frame_len);
out = [];
for i = 1 : frame_num
    iniind = 1 + ( i - 1 ) * frame_shift;
    endind = iniind + frame_len - 1;
    xt = mono_audio( iniind : endind ); % ��֡
    xt_win = wint .* xt; % �Ӵ�
    fft_xt_win = fft( xt_win, nfft ); % FFT
    pf = fft_xt_win( 1 : freq_num ); % ȡ����Ƶ��
    % ��ȫƵ��
    full_pf = [ pf; conj( pf( end - 1 : - 1 : 2 ) ) ];
    % ����Ҷ���任��ʱ��
    inv_full_pf = ifft(full_pf, nfft);
    inv_full_pf_dewin = inv_full_pf( 1 : frame_len ) ./ wint; % ȡ��Ч���ݶ�
    if ( i < frame_num )
        out_xt = inv_full_pf_dewin( 1 : frame_shift );
    else
        out_xt = inv_full_pf_dewin;
    end
    out = [ out; out_xt ];
end
Mo = length( out );
t = ( 1 : M ) / fs;
to = ( 1 : Mo ) / fs;
figure
subplot 211
hold on 
plot(t,mono_audio,'b.');
plot(to,out,'r')
xlim([min(t) max(t)])
legend('ԭʼ�ź�','�Ӵ��ص�����Ļָ��ź�')

subplot 212
hold on 
plot(t,mono_audio,'b.');
plot(to,out,'r')
xlim([4.1 4.15])
legend('ԭʼ�ź�','�Ӵ��ص�����Ļָ��ź�')

