function array_normal = GetArrayNormalIfExists...
    (original_array, element_num, zero_thres)
%如果阵列为线阵或者面阵，则计算阵列法向向量，否则输出零向量

array_direction = GetDirectionIfLinear ...
    (original_array, element_num, zero_thres);

if ( any( array_direction ) ) % 如果阵列为线阵
    array_normal = [array_direction(2),-array_direction(1),0];
else
    array_normal = GetNormalIfPlanar...
        (original_array, element_num, zero_thres);
end








