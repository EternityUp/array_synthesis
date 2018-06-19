#include "fft_transform.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
    bool a = IsPowerOf2(0);
    a = IsPowerOf2(1);
    a = IsPowerOf2(2);
    a = IsPowerOf2(2048);
    a = IsPowerOf2(24);

    size_t b = PowernOf2(0);
    b = PowernOf2(1);
    b = PowernOf2(2);
    b = PowernOf2(10);
    size_t c = NextPower2(0);
    c = NextPower2(1);
    c = NextPower2(2);
    c = NextPower2(123);
    c = NextPower2(257);
    c = NextPower2(2048);


    // 测试傅里叶变换和反变换
    size_t nfft = 2048;
    size_t numFreqBins = nfft / 2 + 1;
    int f0 = 4;  // 信号频率
    int fs = 10; // 采样频率
    float ts = 1.0f / fs; // 采样时间间隔
    float *src = new float[nfft]{0.0f};
    Complexf *dst = new Complexf[numFreqBins];
    float *inv_src = new float[nfft]{0.0f};
    for (size_t i = 0; i < nfft; i++) {
        src[i] = (float) sin(2 * M_PI * f0 * ts * i + M_PI / 4);
    }

    struct timeval start_val, end_val;
    float elapsed_time;
    gettimeofday(&start_val, NULL);
    FftTransform(src, dst, nfft);
    InvFftTransform(dst, inv_src, nfft);
    gettimeofday(&end_val, NULL);
    elapsed_time = (end_val.tv_sec - start_val.tv_sec) * 1000 + 1.f * (end_val.tv_usec - start_val.tv_usec) / 1000.f;

    printf("fft and ifft of %ld points consumes %.4f ms\n", nfft, elapsed_time);
    delete[] src;
    delete[] dst;
    delete[] inv_src;

    return 0;
}