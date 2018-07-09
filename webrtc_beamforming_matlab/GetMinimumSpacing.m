function min_mic_spacing = GetMinimumSpacing(original_array, element_num)
%计算最小的阵元间距
min_mic_spacing = 1e6;
for i = 1 : element_num - 1
    for j = i + 1 : element_num
        spacing = norm(original_array(i,:) - original_array(j,:) );
        min_mic_spacing = min(min_mic_spacing, spacing);
    end
end



