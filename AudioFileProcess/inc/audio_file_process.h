#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_AUDIO_FILE_PROCESS_H_
#define AMY_MODULES_AUDIO_PROCESSING_AUDIO_FILE_PROCESS_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus*/

#ifdef _WIN32
#ifndef AMY_AUDIO_LIB_EXPORTS
#ifdef AMY_AUDIO_EXPORTS 
#define AMY_AUDIO_API _declspec(dllexport)
#else
#define AMY_AUDIO_API _declspec(dllimport)
#endif
#else
#define AMY_AUDIO_API
#endif
#else
#define AMY_AUDIO_API
#endif

#include <stddef.h>
#include <stdint.h>

#define MAX_LIMITS_INT16  32767
#define MIN_LIMITS_INT16  -32768


AMY_AUDIO_API void InputBufferWrite(float *audio_in, float *in_buffer, int chunk_len,
                                    int &read_pos, int &write_pos, bool &rw_wrap, int buffer_len);

AMY_AUDIO_API void InputBufferRead(float *in_buffer, float *in_block, int block_len,
                                   int &read_pos, int &write_pos, bool &rw_wrap, int buffer_len);

AMY_AUDIO_API void MoveReadPositionBackward(int moved_frames, int &read_pos,
                                            const int &write_pos, bool &rw_wrap, int buffer_len);

AMY_AUDIO_API void FloatS16ToFloat(const float *src, size_t size, float *dest);

AMY_AUDIO_API void FloatToFloatS16(const float *src, size_t size, float *dest);

AMY_AUDIO_API void Int16_tToFloat(const int16_t *src, size_t size, float *dest);

AMY_AUDIO_API void FloatToInt16_t(const float *src, size_t size, int16_t *dest);

AMY_AUDIO_API void Deinterleave(const float *interleaved, float **deinterleaved);

AMY_AUDIO_API void Interleave(float **deinterleaved, float *interleaved);

AMY_AUDIO_API void AddFrames(float *buffer, float *block,
                             int buffer_start_ind, int block_start_ind, int len);

AMY_AUDIO_API void OutputSamplesAndUpdateBuffer(float *out_ch_data,
                                                float *buffer, int frame_len, int buffer_len, int initial_dalay);


#ifdef __cplusplus
}
#endif /* __cplusplus*/


#endif











