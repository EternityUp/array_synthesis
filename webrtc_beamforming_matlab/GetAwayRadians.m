function away_radians = GetAwayRadians...
    (kMinAwayRadians, kAwaySlope, min_mic_spacing)
%����������Ԫ��С������س�������Ŀ���������ŷ������С��λ��ֵ
away_radians_compute = kAwaySlope * pi / min_mic_spacing;
away_radians = min( pi, max(kMinAwayRadians,away_radians_compute) );


