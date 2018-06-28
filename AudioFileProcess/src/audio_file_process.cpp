#include "basic_processing_parameters.h"
#include "audio_file_process.h"
#include <math.h>
#include <complex>
#include <assert.h>
#include <string.h>

AMY_AUDIO_API void InputBufferWrite(float *audio_in, float *in_buffer, int chunk_len,
                                    int &read_pos, int &write_pos, bool &rw_wrap, int buffer_len) {
    int reads_available = 0, free_elements = 0, write_elements = 0,
            temp = 0, margin = 0, start_ind = 0, audio_start_ind = 0;
    if (!rw_wrap)
        reads_available = write_pos - read_pos;
    else
        reads_available = buffer_len - read_pos + write_pos;
    free_elements = buffer_len - reads_available;
    write_elements = free_elements < chunk_len ? free_elements : chunk_len;
    temp = write_elements;
    margin = buffer_len - write_pos;
    if (write_elements > margin) {
        start_ind = write_pos;
        for (int i = 0; i < margin; i++) {
            in_buffer[start_ind] = audio_in[i];
            start_ind++;
        }
        write_pos = 0;
        temp -= margin;
        rw_wrap = true;
    }
    start_ind = write_pos;
    audio_start_ind = write_elements - temp;
    for (int i = 0; i < temp; i++) {
        in_buffer[start_ind] = audio_in[audio_start_ind];
        start_ind++;
        audio_start_ind++;
    }
    write_pos += temp;
}

AMY_AUDIO_API void InputBufferRead(float *in_buffer, float *in_block, int block_len,
                                   int &read_pos, int &write_pos, bool &rw_wrap, int buffer_len) {
    int reads_available = 0, read_elements = 0, margin = 0, free_elements = 0;
    int first_block_start_ind = 0, first_block_len = 0;
    int second_block_start_ind = 0, second_block_len = 0;
    int first_audio_start_ind = 0, second_audio_start_ind = 0;
    if (!rw_wrap)
        reads_available = write_pos - read_pos;
    else
        reads_available = buffer_len - read_pos + write_pos;
    read_elements = reads_available < block_len ? reads_available : block_len;
    margin = buffer_len - read_pos;

    if (read_elements > margin) {
        /* write data in two blocks that wrap the buffer */
        first_block_start_ind = read_pos;
        first_block_len = margin;
        second_block_start_ind = 0;
        second_block_len = read_elements - margin;
    } else {
        first_block_start_ind = read_pos;
        first_block_len = read_elements;
        second_block_start_ind = 0;
        second_block_len = 0;
    }

    if (second_block_len > 0) {
        first_audio_start_ind = 0;
        second_audio_start_ind = first_block_len;
        for (int i = 0; i < first_block_len; i++) {
            in_block[first_audio_start_ind] = in_buffer[first_block_start_ind];
            first_audio_start_ind++;
            first_block_start_ind++;
        }
        for (int i = 0; i < second_block_len; i++) {
            in_block[second_audio_start_ind] = in_buffer[second_block_start_ind];
            second_audio_start_ind++;
            second_block_start_ind++;
        }
    } else {
        first_audio_start_ind = 0;
        for (int i = 0; i < first_block_len; i++) {
            in_block[first_audio_start_ind] = in_buffer[first_block_start_ind];
            first_audio_start_ind++;
            first_block_start_ind++;
        }
    }

    free_elements = buffer_len - reads_available;
    if (read_elements > reads_available)
        read_elements = reads_available;
    if (read_elements < -free_elements)
        read_elements = -free_elements;
    read_pos += read_elements;

    if (read_pos > buffer_len) {
        read_pos -= buffer_len;
        rw_wrap = false;
    }

    if (read_pos < 0) {
        read_pos += buffer_len;
        rw_wrap = true;
    }
}

AMY_AUDIO_API void MoveReadPositionBackward(int moved_frames, int &read_pos,
                                            const int &write_pos, bool &rw_wrap, int buffer_len) {
    int reads_available = 0, free_elements = 0;
    if (!rw_wrap)
        reads_available = write_pos - read_pos;
    else
        reads_available = buffer_len - read_pos + write_pos;
    free_elements = buffer_len - reads_available;
    if (moved_frames > reads_available)
        moved_frames = reads_available;
    if (moved_frames < -free_elements)
        moved_frames = -free_elements;

    read_pos += moved_frames;
    if (read_pos > buffer_len) {
        read_pos -= buffer_len;
        rw_wrap = false;
    }

    if (read_pos < 0) {
        read_pos += buffer_len;
        rw_wrap = true;
    }
}


AMY_AUDIO_API void FloatS16ToFloat(const float *src, size_t size, float *dest) {
    static const float kMaxInt16Inverse = 1.f / MAX_LIMITS_INT16;
    static const float kMinInt16Inverse = 1.f / MIN_LIMITS_INT16;

    for (size_t i = 0; i < size; ++i)
        dest[i] = src[i] * (src[i] > 0 ? kMaxInt16Inverse : -kMinInt16Inverse);
}


AMY_AUDIO_API void FloatToFloatS16(const float *src, size_t size, float *dest) {
    for (size_t i = 0; i < size; ++i)
        dest[i] = src[i] * (src[i] > 0 ? MAX_LIMITS_INT16 : - MIN_LIMITS_INT16);
}

AMY_AUDIO_API void Int16_tToFloat(const int16_t *src, size_t size, float *dest) {
    static const float kMaxInt16Inverse = 1.f / MAX_LIMITS_INT16;
    static const float kMinInt16Inverse = 1.f / MIN_LIMITS_INT16;
    for (size_t i = 0; i < size; i++) {
        dest[i] = src[i] * (src[i] > 0 ? kMaxInt16Inverse : -kMinInt16Inverse);
    }
}

AMY_AUDIO_API void FloatToInt16_t(const float *src, size_t size, int16_t *dest) {
    for (size_t i = 0; i < size; ++i) {
        if (src[i] > 0)
            dest[i] = src[i] >= 1 ? MAX_LIMITS_INT16 : (int16_t) (src[i] * MAX_LIMITS_INT16 + 0.5f);
        dest[i] = src[i] <= -1 ? MIN_LIMITS_INT16 : (int16_t) (-src[i] * MIN_LIMITS_INT16 - 0.5f);
    }
}


void Deinterleave(const float *interleaved, float **deinterleaved, size_t channels, size_t len) {
    for (size_t i = 0; i < channels; ++i) {
        float *channel = deinterleaved[i];
        size_t interleaved_idx = i;
        for (size_t j = 0; j < len; ++j) {
            channel[j] = interleaved[interleaved_idx];
            interleaved_idx += channels;
        }
    }
}


void Interleave(float **deinterleaved, float *interleaved, size_t channels, size_t len) {
    for (size_t i = 0; i < channels; ++i) {
        const float *channel = deinterleaved[i];
        size_t interleaved_idx = i;
        for (size_t j = 0; j < len; ++j) {
            interleaved[interleaved_idx] = channel[j];
            interleaved_idx += channels;
        }
    }
}


AMY_AUDIO_API void AddFrames(float *buffer, float *block,
                             int buffer_start_ind, int block_start_ind, int len) {
    for (int j = 0; j < kBlockLen; j++) {
        buffer[buffer_start_ind] += block[block_start_ind];
        buffer_start_ind++;
        block_start_ind++;
    }
}


AMY_AUDIO_API void OutputSamplesAndUpdateBuffer(float *out_ch_data,
                                                float *buffer, int frame_len, int buffer_len, int initial_dalay) {
    // Copy output buffer to output
    int j = 0;
    for (j = 0; j < frame_len; j++)
        out_ch_data[j] = buffer[j];

    // Copy output buffer [frame_len, frame_len + initial_dalay)
    // to output buffer [0, initial_dalay), zero the rest.
    int copy_ind = frame_len;
    for (j = 0; j < initial_dalay; j++) {
        buffer[j] = buffer[copy_ind];
        copy_ind++;
    }
    for (j = initial_dalay; j < buffer_len; j++)
        buffer[j] = 0;
}