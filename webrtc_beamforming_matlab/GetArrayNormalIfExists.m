function array_normal = GetArrayNormalIfExists...
    (original_array, element_num, zero_thres)
%�������Ϊ�������������������з����������������������

array_direction = GetDirectionIfLinear ...
    (original_array, element_num, zero_thres);

if ( any( array_direction ) ) % �������Ϊ����
    array_normal = [array_direction(2),-array_direction(1),0];
else
    array_normal = GetNormalIfPlanar...
        (original_array, element_num, zero_thres);
end








