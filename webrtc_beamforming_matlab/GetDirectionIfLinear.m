function array_direction = GetDirectionIfLinear ...
(original_array, element_num, zero_thres)
%����ж�����Ϊ������������з����������������������

% �����1�ź͵�2����Ԫ�����߷���
first_pair_direction = original_array(2,:) - original_array(1,:);
for i = 2 : element_num
    pair_direction = original_array(i,:) - original_array(i-1,:);
    % �������������Ĳ��
    cross_multiply = cross(first_pair_direction,pair_direction);
    % �����˽���ĵ��
    dot_multiply = dot(cross_multiply,cross_multiply);
    % �����ж�
    if (dot_multiply > zero_thres)
        array_direction = [0.0 0.0 0.0];
    else
        array_direction = first_pair_direction;
    end
end




