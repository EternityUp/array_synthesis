function smoothedPf = smoothPf( pf, winlen, h, M, N )
%��ͨ��Ƶ��ƽ��
% Pf:��ƽ���Ķ�ͨ��Ƶ��
% winlen:ƽ�����ڴ�С
% h:ƽ����������
% M:Ƶ��Ƶ�ʵ���
% N:ͨ����Ŀ
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
