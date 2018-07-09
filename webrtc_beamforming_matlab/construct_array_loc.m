function array_loc = construct_array_loc(M,d,array_type)
%根据阵元数目和阵型构造阵列阵元位置信息
array_loc = zeros(M,3);
if ( array_type == 0 )
    for i = 1 : M
        array_loc(i,1) = ( i - 1 ) * d;
    end
elseif ( array_type == 1 )
    r = d; % 圆阵半径
    theta = 360 / M; % 相邻阵元的圆心角
    array_loc = zeros(M,3);
    for i = 1 : M
        array_loc(i,1) = r * cosd( ( i - 1 ) * theta );
        array_loc(i,2) = r * sind( ( i - 1 ) * theta );
    end
end

