function away_radians = GetAwayRadians...
    (kMinAwayRadians, kAwaySlope, min_mic_spacing)
%根据阵列阵元最小间距和相关常数计算目标和最近干扰分离的最小方位差值
away_radians_compute = kAwaySlope * pi / min_mic_spacing;
away_radians = min( pi, max(kMinAwayRadians,away_radians_compute) );


