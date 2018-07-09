function smoothedPf = smoothPf( pf, winlen, h, M, N )
%多通道频谱平滑
% Pf:待平滑的多通道频谱
% winlen:平滑窗口大小
% h:平滑窗口序列
% M:频谱频率点数
% N:通道数目
smoothedPf = zeros( M, N );
half_len = ( winlen - 1 ) / 2;
for i = 1 : N
    for j = 1 : M
        j_start = max( 1, j - half_len );
        j_end = min( M, j + half_len );
        j_range = j_start : j_end;
        smoothedPf( j, i ) = h(j_range - j_start + 1)' * pf( j_range, i ) ...
            / length(j_range);
       % smoothedPf( j, i ) = mean( pf( j_range, i ) );
    end
end
