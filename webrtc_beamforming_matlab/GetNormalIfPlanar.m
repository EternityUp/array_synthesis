function array_normal = GetNormalIfPlanar...
    (original_array, element_num, zero_thres)
%如果判定阵列为面阵，则输出阵列法向向量，否则输出零向量

% 判定是否为线阵
array_direction = GetDirectionIfLinear ...
(original_array, element_num, zero_thres);
if ( any( array_direction ) ) % 如果为线阵
    array_normal = [0.0 0.0 0.0];
else
    % 判定是否为面阵
    % 计算第1号和第2号阵元的连线方向
    first_pair_direction = original_array(2,:) - original_array(1,:);
    % 计算第2号和第3号阵元的连线方向
    second_pair_direction = original_array(3,:) - original_array(2,:);
    % 依据first_pair_direction和second_pair_direction计算前三个
    % 阵元组成的面阵的法向向量
    normal_direction = cross(first_pair_direction,second_pair_direction);
    % 判断先前计算的法向向量与其他阵元的连线方向是否垂直
    for i = 2 : element_num - 1
        pair_direction = original_array(i+1,:) - original_array(i,:);
        % 计算点积
        dot_multiply = dot(normal_direction,pair_direction);
        if ( dot(dot_multiply,dot_multiply) > zero_thres ) % 如果为立体阵
            array_normal = [0.0 0.0 0.0];
            break;
        else
            array_normal = normal_direction;
        end
    end
end




