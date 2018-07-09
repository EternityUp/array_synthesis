function array_loc = construct_array_loc(M,d,array_type)
%������Ԫ��Ŀ�����͹���������Ԫλ����Ϣ
array_loc = zeros(M,3);
if ( array_type == 0 )
    for i = 1 : M
        array_loc(i,1) = ( i - 1 ) * d;
    end
elseif ( array_type == 1 )
    r = d; % Բ��뾶
    theta = 360 / M; % ������Ԫ��Բ�Ľ�
    array_loc = zeros(M,3);
    for i = 1 : M
        array_loc(i,1) = r * cosd( ( i - 1 ) * theta );
        array_loc(i,2) = r * sind( ( i - 1 ) * theta );
    end
end

