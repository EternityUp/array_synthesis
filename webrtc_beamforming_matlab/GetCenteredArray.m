function centered_array = GetCenteredArray(original_array, element_num)
%�������и���Ԫλ���������ģ��Դ�Ϊ�ο������±궨����
avg_x_y_z = mean(original_array,1);
centered_array = original_array - repmat(avg_x_y_z,element_num,1);







