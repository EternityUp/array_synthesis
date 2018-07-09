function array_normal = GetNormalIfPlanar...
    (original_array, element_num, zero_thres)
%����ж�����Ϊ������������з����������������������

% �ж��Ƿ�Ϊ����
array_direction = GetDirectionIfLinear ...
(original_array, element_num, zero_thres);
if ( any( array_direction ) ) % ���Ϊ����
    array_normal = [0.0 0.0 0.0];
else
    % �ж��Ƿ�Ϊ����
    % �����1�ź͵�2����Ԫ�����߷���
    first_pair_direction = original_array(2,:) - original_array(1,:);
    % �����2�ź͵�3����Ԫ�����߷���
    second_pair_direction = original_array(3,:) - original_array(2,:);
    % ����first_pair_direction��second_pair_direction����ǰ����
    % ��Ԫ��ɵ�����ķ�������
    normal_direction = cross(first_pair_direction,second_pair_direction);
    % �ж���ǰ����ķ���������������Ԫ�����߷����Ƿ�ֱ
    for i = 2 : element_num - 1
        pair_direction = original_array(i+1,:) - original_array(i,:);
        % ������
        dot_multiply = dot(normal_direction,pair_direction);
        if ( dot(dot_multiply,dot_multiply) > zero_thres ) % ���Ϊ������
            array_normal = [0.0 0.0 0.0];
            break;
        else
            array_normal = normal_direction;
        end
    end
end




