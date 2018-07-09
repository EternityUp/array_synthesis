function [input_block_,read_pos,write_pos,rw_wrap] = ...
    InputBufferRead(input_buffer_,input_block_,block_size_,...
    read_pos,write_pos,rw_wrap,element_count)
%将数据从input_buffer_写入input_block_中，同时更新读指针、写指针的位置以及读写指针的
%相对位置状态
if (rw_wrap == 0 )
    reads_available = write_pos - read_pos;
else
    reads_available = element_count - read_pos + write_pos;
end
read_elements = min(reads_available, block_size_);
margin = element_count - read_pos;

if ( read_elements > margin )
    % write data in two blocks that wrap the buffer
    block_start_index1 = read_pos + 1;
    block_len1 = margin;
    block_start_index2 = 1;
    block_len2 = read_elements - margin;
else
    block_start_index1 = read_pos + 1;
    block_len1 = read_elements;
    block_start_index2 = 0;
    block_len2 = 0;
end

if( block_len2 > 0 )
    source_index_range1 = block_start_index1 : ... 
        (block_start_index1 + block_len1 - 1);
    source_index_range2 = block_start_index2 : ... 
        (block_start_index2 + block_len2 - 1);
    dst_index_range1 = 1 : block_len1;
    dst_index_range2 = ( block_len1 + 1 ) : ( block_len1 + block_len2 );
    input_block_(dst_index_range1,:) = ... 
        input_buffer_(source_index_range1,:);
    input_block_(dst_index_range2,:) = ... 
        input_buffer_(source_index_range2,:);
else
    source_index_range = block_start_index1 : ... 
        (block_start_index1 + block_len1 - 1);
    dst_index_range = 1 : block_len1;
    input_block_(dst_index_range,:) = input_buffer_(source_index_range,:);
end

% 更新度指针位置
free_elements = element_count - reads_available;
readable_elements = reads_available;
if ( read_elements > readable_elements )
    read_elements = readable_elements;
end
if ( read_elements < -free_elements )
    read_elements = -free_elements;
end
read_pos = read_pos + read_elements;
if ( read_pos > element_count )
    read_pos = read_pos - element_count;
    rw_wrap = 0;
end

if( read_pos < 0 )
    read_pos = read_pos + element_count;
    rw_wrap = 1;
end


