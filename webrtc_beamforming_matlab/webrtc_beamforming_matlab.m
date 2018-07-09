clear all
close all
clc
tic
%% 读取文件
% [audio,fs] = audioread('stereo.wav');
% [audio,fs] = audioread('cos_stereo.wav');
% [audio,fs] = audioread('stereo_linear.wav');
% [audio,fs] = audioread('stereo_circular.wav');
[audio,fs] = audioread('xmos_record2_chs4.wav');
MM = size(audio,2);
%% 常量
kChunksPerSecond = 100; % 每秒包含的语音块数目
kChunkSizeMs = 1000 / kChunksPerSecond; % 每帧语音时长


num_input_channels = MM; % 输入通道数目
num_postfilter_channels  = MM; % 后置滤波通道数目
% original_array_geometry = [ ...
%     0.00 0.00 0.00;...
%     0.05 0.00 0.00;...
%     ]; % 原始阵列阵元坐标
% original_array_geometry = construct_array_loc...
%     (MM,0.05,1); % 构造阵列阵元位置

load  array_loc_xmos.mat
original_array_geometry = array_loc_xmos;



element_num = size(original_array_geometry,1); % 阵元数目
kMaxDotProduct = 1e-6; % 零值判定阈值
% 计算重新标定后的阵列阵元位置坐标
centered_array = GetCenteredArray(original_array_geometry,element_num);
% 判断阵型
line_direction = GetDirectionIfLinear ...
    (original_array_geometry, element_num, kMaxDotProduct);
planar_normal = GetNormalIfPlanar...
    (original_array_geometry, element_num, kMaxDotProduct);
array_normal = GetArrayNormalIfExists...
    (original_array_geometry, element_num, kMaxDotProduct);
% 计算最小的阵元间距
min_mic_spacing = GetMinimumSpacing(original_array_geometry, element_num);
de2rad = pi / 180;
rad2de = 180 / pi;
target_angle_radians_ = 60 * de2rad; % 目标方位
kMinAwayRadians = 0.2; % 认为目标和最近干扰方位分离的最小角度差(rad)
kAwaySlope = 0.008; % 影响目标和最近干扰方位分离的最小角度差的常数因子
% 计算目标和干扰分离的最小角度差(rad)
away_radians = GetAwayRadians...
    (kMinAwayRadians, kAwaySlope, min_mic_spacing);
% 计算窗函数
kKbdAlpha = 1.5; % Kaiser Bessel Derived window参数
kFftSize = 256;
kNumFreqBins = kFftSize / 2 + 1;
wint = KaiserBesselDerived(kKbdAlpha, kFftSize);

%% 其他相关常量值
kSpeedOfSoundMeterSeconds = 343; % 空气中声速
% 计算干扰协方差矩阵时的加权因子
% Rpsi = Rpsi_angled * kBalance + Rpsi_uniform * (1 - kBalance)
kBalance = 0.95;
% 计算掩蔽值的平滑因子
kMaskTimeSmoothAlpha = 0.2; % 时域
kMaskFrequencySmoothAlpha = 0.6; % 频域
% 计算中频段的平均掩蔽值，与kMaskQuantile有关
kLowMeanStartHz = 200;
kLowMeanEndHz = 400;
kCutOffConstant = 0.9999;
kMaskQuantile = 0.7;
kMaskTargetThreshold = 0.01;
kHoldTargetSeconds = 0.25;
kCompensationGain = 2.0;

%% bf相关初始化操作
chunk_length_ = round(fs / kChunksPerSecond); % 每帧数据长度
high_pass_postfilter_mask_ = 1.0;
is_target_present_ = 0;
hold_target_blocks_ = round(kHoldTargetSeconds * 2 * fs / kFftSize);
interference_blocks_count_ = hold_target_blocks_;
wave_number_step = ...
    (2.0 * pi * fs) / (kFftSize * kSpeedOfSoundMeterSeconds);

time_smooth_mask_ = zeros(kNumFreqBins,1);
final_mask_ = zeros(1,kNumFreqBins);
wave_numbers_ = zeros(1,kNumFreqBins);
for i = 1 : kNumFreqBins
    time_smooth_mask_(i) = 1.0;
    final_mask_(i) = 1.0;
    wave_numbers_(i) = ( i - 1 ) * wave_number_step;
end
%% 计算低频段频率下标范围
low_mean_start_bin_ = round( kLowMeanStartHz * ...
    kFftSize / fs) + 1;
low_mean_end_bin_ = round( kLowMeanEndHz * ...
    kFftSize / fs) + 1;
%% 计算阵元之间的间距
D = zeros(element_num);
for d1 = 1 : element_num
    for d2 = 1 : element_num
        D(d1,d2) = norm( centered_array(d1, :) - ...
            centered_array(d2, :) );
    end
end
%% 计算频域散射协方差矩阵
uniform_cov_mat_ = zeros(kNumFreqBins,element_num,element_num);
for i = 1 : kNumFreqBins
    if( wave_numbers_(i) > 0 )
        uniform_cov_mat_(i,:,:) = besselj(0, wave_numbers_(i) * D );
    else
        uniform_cov_mat_(i,:,:) = eye(element_num);
    end
    normalization_factor = uniform_cov_mat_(i,1,1);
    uniform_cov_mat_(i,:,:) = uniform_cov_mat_(i,:,:) ...
        / normalization_factor * ( 1 - kBalance );
end
%% 计算混叠频率和高频段频率下标范围
kAliasingFreqHz = kSpeedOfSoundMeterSeconds / ...
    (min_mic_spacing * (1.0 + abs(cos(target_angle_radians_))));
kHighMeanStartHz = min(0.5 *  kAliasingFreqHz, fs / 2.0);
kHighMeanEndHz = min(0.75 *  kAliasingFreqHz, fs / 2.0);
high_mean_start_bin_ = round(kHighMeanStartHz * kFftSize / fs) + 1;
high_mean_end_bin_ = round(kHighMeanEndHz * kFftSize / fs) + 1;
%% 计算干扰方向(顺时针和逆时针)
diff_rad = target_angle_radians_ - away_radians;
target_direction = [cos(target_angle_radians_) ...
    sin(target_angle_radians_) 0];
clockwise_interf_direction = [cos(diff_rad) sin(diff_rad) 0];
if ( any(array_normal) || ...
        dot(array_normal, target_direction) * ...
        dot(array_normal, clockwise_interf_direction) >= 0.0)
    interf_angles_radians_(1) = target_angle_radians_ - away_radians;
else
    interf_angles_radians_(1) = target_angle_radians_ - away_radians + pi;
end
diff_rad_counterclock = target_angle_radians_ + away_radians;
counterclock_interf_direction = [cos(diff_rad_counterclock) ...
    sin(diff_rad_counterclock) 0];
if ( any(array_normal) || ...
        dot(array_normal, target_direction) * ...
        dot(array_normal, clockwise_interf_direction) >= 0.0)
    interf_angles_radians_(2) = target_angle_radians_ + away_radians;
else
    interf_angles_radians_(2) = target_angle_radians_ + away_radians - pi;
end
Li = length(interf_angles_radians_);
%% 计算时延导向矢量
delay_sum_masks_ = zeros(kNumFreqBins,element_num);
% 计算声程差(只考虑x-y平面)
target_direction2D = target_direction(1:2);
sound_path_difference = target_direction2D * centered_array(:,1:2)';
for i = 1 : kNumFreqBins
    freq_in_hertz = (i - 1) / kFftSize * fs;
    phase_shift = -2.0 * pi * sound_path_difference * freq_in_hertz ...
        / kSpeedOfSoundMeterSeconds;
    delay_sum_masks_(i,:) = exp( 1j * phase_shift );
    norm_factor = norm(delay_sum_masks_(i,:),2);
    delay_sum_masks_(i,:) = delay_sum_masks_(i,:) / norm_factor;
end
%% 计算目标方向协方差矩阵
target_cov_mats_ = zeros(kNumFreqBins, element_num, element_num);
for i = 1 : kNumFreqBins
    target_cov_mats_(i,:,:) = delay_sum_masks_(i,:).' * ...
        conj( delay_sum_masks_(i, :) );
end
%% 计算干扰方向协方差矩阵
interf_cov_mats_ = zeros(kNumFreqBins, Li, element_num, element_num);
% 计算干扰方向声程差
interferer_path_difference = zeros(Li,element_num);
for i = 1 : Li
    interferer_doa = interf_angles_radians_(i);
    interferer_direction2D = [cos(interferer_doa) ...
        sin(interferer_doa)];
    interferer_path_difference(i,:) = interferer_direction2D * ...
        centered_array(:,1:2)';
end
for i = 1 : kNumFreqBins
    freq_in_hertz = ( i - 1 ) / kFftSize * fs;
    for j = 1 : Li
        interferer_phase_shift = -2.0 * pi * interferer_path_difference(j,:)...
            * freq_in_hertz / kSpeedOfSoundMeterSeconds;
        interferer_delay_sum_masks_ = ...
            exp( 1j * interferer_phase_shift );
        interferer_delay_sum_masks_ = interferer_delay_sum_masks_ ...
            / norm(interferer_delay_sum_masks_);
        interferer_cov_matrix = interferer_delay_sum_masks_.' * ...
            conj(interferer_delay_sum_masks_);
        normalization_factor = interferer_cov_matrix(1,1);
        interferer_cov_matrix = interferer_cov_matrix ...
            / normalization_factor * kBalance;
        interf_cov_mats_(i,j,:,:) = interferer_cov_matrix + ...
            squeeze(uniform_cov_mat_(i,:,:));
    end
end
%% 计算归一化波束方位能量(目标和干扰方位)
rxiws_ = zeros(kNumFreqBins,1);
rpsiws_ = zeros(kNumFreqBins,2);
for i = 1 : kNumFreqBins
    rxiws_(i) = conj(delay_sum_masks_(i,:)) * squeeze(target_cov_mats_(i,:,:)) * ...
        delay_sum_masks_(i,:).';
    for j = 1 : Li
        rpsiws_(i,j) = conj(delay_sum_masks_(i,:)) * squeeze(interf_cov_mats_(i,j,:,:)) * ...
            delay_sum_masks_(i,:).';
    end
end
%% 分帧处理
L = length(audio);
frame_offset_ = 0;
block_size_ = 256;
shift_amount_ = 128;
initial_delay_ =  block_size_ - ( chunk_length_ - shift_amount_ );
element_count = chunk_length_ + initial_delay_;
input_buffer_ = zeros( element_count, MM );
output_buffer_ = zeros( element_count, MM );
input_block_ = zeros( block_size_, MM );
output_block_ = zeros( block_size_, MM );
read_pos = 160;
write_pos = 0;
rw_wrap = 1; % 0:SAME_WRAP   1:DIFF_WRAP
CASE = 0;
frame_num = floor( L / chunk_length_ );
new_mask_ = zeros(kNumFreqBins,1);
audio_out = [];
for i = 1 : frame_num
    prev_fft_block = [];
    first_frame_in_block_list = [];
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind, : );
    [input_buffer_,read_pos,write_pos,rw_wrap] = ...
        InputBufferWrite(buf,input_buffer_,chunk_length_,...
        read_pos,write_pos,rw_wrap,element_count);
    first_frame_in_block = frame_offset_;
    %% AnalyzeChunk
    while ( first_frame_in_block < chunk_length_ )
        [input_block_,read_pos,write_pos,rw_wrap] = ...
            InputBufferRead(input_buffer_,input_block_,block_size_,...
            read_pos,write_pos,rw_wrap,element_count);
        moved_frames = -block_size_ + shift_amount_;
        [read_pos,write_pos,rw_wrap] = ...
            InputBufferMoveReadPositionBackward(moved_frames,read_pos,...
            write_pos,rw_wrap,element_count);
        % 加窗
        input_block_ = input_block_ .* repmat(wint,1,element_num);
        % 频域变换
        fft_block = fft(input_block_,kFftSize,1);
        fft_block = fft_block(1:kNumFreqBins,:);
        prev_fft_block = [prev_fft_block;fft_block];
        % 频域处理，计算后置滤波器掩蔽值
        for k = low_mean_start_bin_ : high_mean_end_bin_
            eig_m_ = fft_block(k,:);
            eig_m_norm_factor = norm(eig_m_);
            if ( eig_m_norm_factor ~= 0 )
                eig_m_ = eig_m_ / eig_m_norm_factor;
            end
            rxim = conj(eig_m_) * squeeze(target_cov_mats_(k,:,:)) * eig_m_.';
            ratio_rxiw_rxim = 0.0;
            if ( rxim > 0.0 )
                ratio_rxiw_rxim = rxiws_(k) / rxim;
            end
            rmw = abs( dot( delay_sum_masks_(k,:), eig_m_ ) ) .^ 2;
            rmw_r = real(rmw);
            % 计算后置滤波器掩蔽值
            % 某频点、第一个干扰方向的干扰协方差矩阵
            interf_cov_mats_fd = squeeze(interf_cov_mats_(k,1,:,:));
            new_mask_(k) = CalculatePostfilterMask(interf_cov_mats_fd,...
                rpsiws_(k,1),ratio_rxiw_rxim,rmw_r,eig_m_,kCutOffConstant);
            for m = 2 : Li
                interf_cov_mats_fd = squeeze(interf_cov_mats_(k,m,:,:));
                tmp_mask = CalculatePostfilterMask(interf_cov_mats_fd,...
                    rpsiws_(k,m),ratio_rxiw_rxim,rmw_r,eig_m_,kCutOffConstant);
            end
            if ( tmp_mask < new_mask_(k) )
                new_mask_(k) = tmp_mask;
            end
            % 掩蔽值平滑
            time_smooth_mask_(k) = kMaskTimeSmoothAlpha * new_mask_(k) + ...
                (1 - kMaskTimeSmoothAlpha) * time_smooth_mask_(k);
        end  % loop mask
        % 判断目标是否存在
        quantile = floor((high_mean_end_bin_ - low_mean_start_bin_) * kMaskQuantile ...
            + low_mean_start_bin_);
        sorted_new_mask_ = ... 
            sort(new_mask_(low_mean_start_bin_ : high_mean_end_bin_),'ascend');
        if ( sorted_new_mask_( quantile - low_mean_start_bin_ + 1) > kMaskTargetThreshold )
            is_target_present_ = 1;
            interference_blocks_count_ = 0;
        else
            interference_blocks_count_ = interference_blocks_count_ + 1;
            is_target_present_ = interference_blocks_count_ < hold_target_blocks_;
        end
        %% 掩蔽值修正
        % 低频mask值修正
        low_freq_range = low_mean_start_bin_ : low_mean_end_bin_;
        low_frequency_mask = mean(time_smooth_mask_(low_freq_range));
        time_smooth_mask_( 1 : ( low_mean_start_bin_ - 1 ) ) = low_frequency_mask;
        % 高频mask值修正
        high_freq_range = high_mean_start_bin_ : high_mean_end_bin_;
        high_pass_postfilter_mask_ = mean(time_smooth_mask_(high_freq_range));
        time_smooth_mask_( ( high_mean_end_bin_ + 1 ) : end ) = high_pass_postfilter_mask_;
        
        % 频域平滑处理 ---> final_mask
        final_mask_ = time_smooth_mask_;
        for k1 = low_mean_start_bin_ : kNumFreqBins
            final_mask_(k1) = kMaskFrequencySmoothAlpha * final_mask_(k1) + ...
            (1 - kMaskFrequencySmoothAlpha) * final_mask_(k1-1);
        end
        for k2 = ( high_mean_end_bin_ ) + 1 : -1 : 2
            final_mask_(k2-1) = kMaskFrequencySmoothAlpha * final_mask_(k2-1) + ...
            (1 - kMaskFrequencySmoothAlpha) * final_mask_(k2);
        end
        first_frame_in_block_list = [first_frame_in_block_list ...
            first_frame_in_block];
        first_frame_in_block = first_frame_in_block + shift_amount_;
    end
    %% PostFilter
    [M1,N1] = size(prev_fft_block);
    post_fft_block = zeros(M1,N1);
    num_sections = M1 / kNumFreqBins;
    for n1 = 1 : num_sections
        fred_ind1 = 1 + ( n1 - 1 ) * kNumFreqBins;
        fred_ind2 = fred_ind1 + kNumFreqBins - 1;
        rangef = fred_ind1 : fred_ind2;
        input_fft_block = prev_fft_block(rangef,:);
        output_fft_block = zeros(kNumFreqBins,element_num);
        for n2 = 1 : kNumFreqBins
            % 对原始频谱进行掩膜处理
            output_fft_block(n2,:) = kCompensationGain * final_mask_(n2)...
                * input_fft_block(n2,:);
        end
        % 补全频谱
        full_output_fft_block = [output_fft_block;... 
            conj(output_fft_block(end-1:-1:2,:))];
        % 逆傅里叶变换回时域
        output_block_ = ifft(full_output_fft_block,[],1);
        % 加窗
        output_block_ = output_block_ .* repmat(wint,1,element_num);
        % AddFrames
        first_frame_in_blocki = first_frame_in_block_list(n1);
        result_ind_range = (first_frame_in_blocki + 1) : ... 
            (first_frame_in_blocki + block_size_);
        b_ind_range = 1 : block_size_;
        a_ind_range = (first_frame_in_blocki + 1) : ... 
            (first_frame_in_blocki + block_size_);
        output_buffer_(result_ind_range, :) =  ...
            output_buffer_(a_ind_range, :) + ...
            output_block_(b_ind_range,:);
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
    
    % 计算新的起始帧偏移
    frame_offset_ = first_frame_in_block - chunk_length_;
    audio_out = [audio_out;out_buf];
end
norm_factor = max(max(abs(audio_out)));
audio_out = audio_out / norm_factor;
audiowrite('stereo2_out.wav',audio_out,fs);
toc
