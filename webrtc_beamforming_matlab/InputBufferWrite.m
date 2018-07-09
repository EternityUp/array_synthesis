function [input_buffer_,read_pos,write_pos,rw_wrap] = ...
    InputBufferWrite(buf,input_buffer_,chunk_length_,...
    read_pos,write_pos,rw_wrap,element_count)
%将数据从buf写入input_buffer中，同时更新读指针、写指针的位置以及读写指针的
%相对位置状态

if ( rw_wrap == 0 )
    reads_available = write_pos - read_pos;
else
    reads_available = element_count - read_pos + write_pos;
end
free_elements = element_count - reads_available;
write_elements = min(free_elements, chunk_length_);
n = write_elements;
margin = element_count - write_pos;
if ( write_elements > margin )
    ind_range = ( write_pos + 1 ) : ( write_pos + margin );
    input_buffer_( ind_range, : ) = buf( 1 : margin, : );
    write_pos = 0;
    n = n - margin;
    rw_wrap = 1;
end
ind_range = ( write_pos + 1 ) : ( write_pos + n );
ind_range2 = ( write_elements - n + 1 ) : write_elements;
input_buffer_( ind_range, : ) = buf( ind_range2, : );
write_pos  = write_pos + n;


