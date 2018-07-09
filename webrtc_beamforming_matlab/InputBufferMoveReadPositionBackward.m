function [read_pos,write_pos,rw_wrap] = ...
    InputBufferMoveReadPositionBackward(moved_frames,read_pos,... 
    write_pos,rw_wrap,element_count)
%将读指针向后移动适当的位置
if (rw_wrap == 0 )
    reads_available = write_pos - read_pos;
else
    reads_available = element_count - read_pos + write_pos;
end
free_elements = element_count - reads_available;
readable_elements = reads_available;
if ( moved_frames > readable_elements )
    moved_frames = readable_elements;
end

if ( moved_frames < -free_elements )
    moved_frames = -free_elements;
end

read_pos = read_pos + moved_frames;
if ( read_pos > element_count )
    read_pos = read_pos - element_count;
    rw_wrap = 0;
end

if( read_pos < 0 )
    read_pos = read_pos + element_count;
    rw_wrap = 1;
end








