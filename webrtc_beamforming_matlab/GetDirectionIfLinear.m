function array_direction = GetDirectionIfLinear ...
(original_array, element_num, zero_thres)
%如果判定阵列为线阵，则输出阵列方向向量，否则输出零向量

% 计算第1号和第2号阵元的连线方向
first_pair_direction = original_array(2,:) - original_array(1,:);
for i = 2 : element_num
    pair_direction = original_array(i,:) - original_array(i-1,:);
    % 计算两个向量的叉积
    cross_multiply = cross(first_pair_direction,pair_direction);
    % 计算叉乘结果的点积
    dot_multiply = dot(cross_multiply,cross_multiply);
    % 线阵判定
    if (dot_multiply > zero_thres)
        array_direction = [0.0 0.0 0.0];
    else
        array_direction = first_pair_direction;
    end
end




