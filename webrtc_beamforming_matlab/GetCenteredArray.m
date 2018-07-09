function centered_array = GetCenteredArray(original_array, element_num)
%计算阵列各阵元位置坐标中心，以此为参考，重新标定坐标
avg_x_y_z = mean(original_array,1);
centered_array = original_array - repmat(avg_x_y_z,element_num,1);







