#include "basic_signal_processing.h"
#include <stdlib.h>
#define LEN 6
int main(int argc, char * argv[])
{
	int length = 16;
	float *Hanning_win = new float[length];
	float *Kaiser_win = new float[length];
	float alpha = 1.5f;

	Hanning(length, Hanning_win);
	KaiserBesselDerived(alpha, length, Kaiser_win);

	delete[] Hanning_win;
	delete[] Kaiser_win;

	int16_t vector1[LEN] = { 1,2,3,2,1,4 };
	int16_t vector2[LEN] = { 0,2,1,1,1,3 };
	int scaling = 1;
	int32_t dpws = DotProductWithScale(vector1, vector2, LEN, scaling);
	
	float src[10] = { 1,2,3,4,5,6,7,8,9,10 };
	float dst[10] = { 0 };
	float win[3] = { 0 };
	Hanning(3, win);
	Smooth1Ddata(src, dst, 10, 3, win);

	float norm = VectorNorm(src, 10);
	ScaleVector(src, dst, 10, 1.f / norm);
	
	size_t gcd = ComputeGcd(100, 24);

	/* 构造数据流 */
	float audio[3200] = { 0 };
	for (int i = 0; i < 3200; i++)
		audio[i] = i * 1.0f;
	int chunk_length = 160;
	int frame_num = 3200 / chunk_length;
	int block_len = 256;
	int shift_amount = 128;
	int initial_delay = block_len - ComputeGcd(chunk_length, shift_amount);
	int frame_offset = 0;
	int buffer_len = chunk_length + initial_delay;
	int read_pos = 160;
	int write_pos = 0;
	bool rw_wrap = true;
	int start_ind = 0;
	float buf[160] = { 0 };
	float in_buffer[384] = { 0 };
	float in_block[256] = { 0 };
	int first_frame_in_block = 0;
	int moved_frames = 0;
	// 模拟数据流读取
	for (int i = 0; i < frame_num; i++)
	{
		start_ind = i * chunk_length;
		for (int j = 0; j < chunk_length; j++)
		{
			buf[j] = audio[start_ind];
			start_ind++;
		}
		InputBufferWrite(buf, in_buffer, chunk_length,
			read_pos, write_pos, rw_wrap, buffer_len);
		first_frame_in_block = frame_offset;
		while (first_frame_in_block < chunk_length)
		{
			InputBufferRead(in_buffer, in_block, block_len,
				read_pos, write_pos, rw_wrap, buffer_len);
			moved_frames = -block_len + shift_amount;
			MoveReadPositionBackward(moved_frames, read_pos,
				write_pos, rw_wrap, buffer_len);
			first_frame_in_block += shift_amount;
		}
		frame_offset = first_frame_in_block - chunk_length;
	}



	system("pause");
	return 0;
}