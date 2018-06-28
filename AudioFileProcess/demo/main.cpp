#include <iostream>
#include <string.h>
#include "basic_processing_parameters.h"
#include "audio_file_process.h"
#include <stdio.h>

int main(int argc, char **argv) {
    using std::cout;
    using std::endl;
    cout << "sizeof(float) = " << sizeof(float) << endl;
    cout << "sizeof(double) = " << sizeof(double) << endl;
    cout << "sizeof(int) = " << sizeof(int) << endl;
    cout << "sizeof(short) = " << sizeof(short) << endl;
    cout << "sizeof(long) = " << sizeof(long) << endl;

    FILE *in_fp, *out_fp;
    int16_t frame_buf_int16[kNumChannels * kFrameLen] = {0};
    float frame_buf_float[kNumChannels * kFrameLen] = {0};
    int16_t out_frame_buf_int16[kNumChannels * kFrameLen] = {0};
    float out_frame_buf_float[kNumChannels * kFrameLen] = {0};

    float **in_multi_chs_buf = new float *[kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
        in_multi_chs_buf[i] = new float[kFrameLen]{0};

    float **out_multi_chs_buf = new float *[kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
        out_multi_chs_buf[i] = new float[kFrameLen]{0};

    int read_size = 0, write_size = 0;

    in_fp = fopen(argv[1], "rb");
    out_fp = fopen(argv[2], "wb");
    if ((NULL == in_fp) || (NULL == out_fp)) {
        cout << "fopen error" << endl;
        return -1;
    }

    unsigned int frame_num = 0;
    while (!feof(in_fp)) {
        frame_num++;
        read_size = fread(frame_buf_int16, sizeof(int16_t), kNumChannels * kFrameLen, in_fp);
        if ((read_size < 0) || (read_size % kNumChannels != 0)) {
            cout << "fread error!  " << "return value is %lu" << read_size << endl;
            return -1;
        }
        Int16_tToFloat(frame_buf_int16, read_size, frame_buf_float);
        Deinterleave(frame_buf_float, in_multi_chs_buf, kNumChannels, kFrameLen);

        for (int j = 0; j < kNumChannels; j++)
            memcpy(out_multi_chs_buf[j], in_multi_chs_buf[j], sizeof(float) * read_size / kNumChannels);

        Interleave(out_multi_chs_buf, out_frame_buf_float, kNumChannels, kFrameLen);
        FloatToInt16_t(out_frame_buf_float, read_size, out_frame_buf_int16);
        write_size = fwrite(out_frame_buf_int16, sizeof(int16_t), read_size, out_fp);
        if (write_size < 0) {
            cout << "fwrite error!  " << "return value is %lu" << write_size << endl;
            return -1;
        }
        cout << "frame_num = " << frame_num << "  read_size = " << read_size
             << "  write_size = " << write_size << endl;
    }

    fclose(in_fp);
    fclose(out_fp);

    for (int i = 0; i < kNumChannels; i++)
        delete[] in_multi_chs_buf[i];
    delete[] in_multi_chs_buf;


    for (int i = 0; i < kNumChannels; i++)
        delete[] out_multi_chs_buf[i];
    delete[] out_multi_chs_buf;


    return 0;
}


