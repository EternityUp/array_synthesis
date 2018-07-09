%模拟数据流读取
clear all
close all
clc
chunk_length_ = 160;
frame_offset_ = 0;
block_size_ = 256;
shift_amount_ = 128;
initial_delay_ =  block_size_ - gcd( chunk_length_, shift_amount_ );
element_count = chunk_length_ + initial_delay_;
input_buffer_ = zeros( element_count, 1 );
output_buffer_ = zeros( element_count, 1 );
input_block_ = zeros( block_size_, 1 );
output_block_ = zeros( block_size_, 1 );
read_pos = 160;
write_pos = 0;
rw_wrap = 1; % 0:SAME_WRAP   1:DIFF_WRAP



audio = ( 0 : 3199 )';
L = length(audio);
frame_num = floor( L / chunk_length_ );
out_buf_course = [];
for i = 1 : frame_num
    ini_ind = 1 + ( i - 1 ) * chunk_length_;
    end_ind = ini_ind + chunk_length_ - 1;
    audio_ind = ini_ind : end_ind;
    buf = audio( audio_ind);
    [input_buffer_,read_pos,write_pos,rw_wrap] = ...
        InputBufferWrite(buf,input_buffer_,chunk_length_,...
        read_pos,write_pos,rw_wrap,element_count);
    first_frame_in_block = frame_offset_;
    while ( first_frame_in_block < chunk_length_ )
        [input_block_,read_pos,write_pos,rw_wrap] = ...
            InputBufferRead(input_buffer_,input_block_,block_size_,...
            read_pos,write_pos,rw_wrap,element_count);
        moved_frames = -block_size_ + shift_amount_;
        [read_pos,write_pos,rw_wrap] = ...
            InputBufferMoveReadPositionBackward(moved_frames,read_pos,...
            write_pos,rw_wrap,element_count);
        
        output_block_ = input_block_;
        % AddFrames
        first_frame_in_blocki = first_frame_in_block;
        result_ind_range = (first_frame_in_blocki + 1) : ...
            (first_frame_in_blocki + block_size_);
        b_ind_range = 1 : block_size_;
        a_ind_range = (first_frame_in_blocki + 1) : ...
            (first_frame_in_blocki + block_size_);
        output_buffer_(result_ind_range, :) =  ...
            output_buffer_(a_ind_range, :) + ...
            output_block_(b_ind_range,:);
        first_frame_in_block = first_frame_in_block + shift_amount_;
    end
    
    % Copy output buffer to output
    % CopyFrames
    out_buf = output_buffer_(1:chunk_length_,:);
    % MoveFrames ----> output_buffer
    dst_ind_range = 1 : initial_delay_;
    src_ind_range = ( chunk_length_ + 1 ) : ...
        ( chunk_length_ + initial_delay_ );
    output_buffer_(dst_ind_range,:) = output_buffer_(src_ind_range,:);
    % ZeroOut ----> output_buffer
    src_ind_range = ( initial_delay_ + 1 ) : ...
        ( chunk_length_ + initial_delay_ );
    output_buffer_(src_ind_range,:) = 0;
    frame_offset_ = first_frame_in_block - chunk_length_;
    out_buf_course = [out_buf_course;out_buf];
end


figure
plot(audio)
hold on
plot(out_buf_course(225:end),'r')
grid on
axis tight
legend('模拟输入序列','未经算法处理的输出序列')