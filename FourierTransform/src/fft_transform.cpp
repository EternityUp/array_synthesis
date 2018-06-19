#include "fft_transform.h"
#include <assert.h>
#include <math.h>
#include <string.h>

AMY_AUDIO_API bool IsPowerOf2(size_t n) {
    int ret;
    if (n == 0)
        return false;
    if (n == 1)
        return true;
    while (n >= 2) {
        ret = n % 2;
        if (ret == 1)
            return false;
        n /= 2;
    }
    return true;
}

AMY_AUDIO_API size_t PowernOf2(size_t n) {
    if (n == 0)
        return 1;
    size_t tmp = 1;
    while (n-- > 0)
        tmp *= 2;
    return tmp;
}

AMY_AUDIO_API size_t NextPower2(size_t n) {
    if (n <= 1)
        return 0;
    size_t tmp = 1;
    while (n > PowernOf2(tmp))
        tmp++;
    return tmp;

}

AMY_AUDIO_API void FftTransform(const float *const in_td_data, Complexf *out_fd_data, size_t len) {
    size_t powern = NextPower2(len);
    size_t nfft = PowernOf2(powern);
    assert(len == nfft && len >= 2);

    size_t half_nfft = nfft / 2;
    size_t ip_length = 2 + int(sqrtf(nfft * 1.0f));
    size_t w_length = half_nfft + 1;

    size_t *ip = new size_t[ip_length]{0};
    float *w = new float[w_length]{0.0f};
    float *tmp = new float[nfft]{0.0f};
    memcpy(tmp, in_td_data, nfft * sizeof(float));

    rdft(nfft, 1, tmp, ip, w);

    out_fd_data[half_nfft] = Complexf(tmp[1], 0.0f);
    out_fd_data[0] = Complexf(tmp[0], 0.0f);

    for (size_t i = 1; i < half_nfft; i++)
        out_fd_data[i] = Complexf(tmp[2 * i], -tmp[2 * i + 1]);


    delete[] ip;
    delete[] w;
    delete[] tmp;
}
// 
AMY_AUDIO_API void InvFftTransform(const Complexf *const in_fd_data, float *out_td_data, size_t len) {
    size_t powern = NextPower2(len);
    size_t nfft = PowernOf2(powern);
    assert(len == nfft && len >= 2);

    size_t half_nfft = nfft / 2;
    size_t ip_length = 2 + int(sqrtf(nfft * 1.0f));
    size_t w_length = nfft / 2 + 1;
    size_t *ip = new size_t[ip_length]{0};
    float *w = new float[w_length]{0.0f};

    out_td_data[1] = in_fd_data[half_nfft].real() / half_nfft;

    for (size_t i = 0; i < half_nfft; i++)
        out_td_data[2 * i] = in_fd_data[i].real() / half_nfft;
    for (size_t i = 1; i < half_nfft; i++)
        out_td_data[2 * i + 1] = -in_fd_data[i].imag() / half_nfft;


    rdft(nfft, -1, out_td_data, ip, w);

    delete[]ip;
    delete[]w;
}


AMY_AUDIO_API void ConstructComplexSpectrum(Complexf *complex_spectrum, float *power_spectrum,
                                            float *phase_spectrum, int nfreq_bins) {
    float real_part = 0.0f, imag_part = 0.0f;
    for (int j = 0; j < nfreq_bins; j++) {
        real_part = sqrtf(power_spectrum[j]) * cosf(phase_spectrum[j]);
        imag_part = sqrtf(power_spectrum[j]) * sinf(phase_spectrum[j]);
        complex_spectrum[j] = Complexf(real_part, imag_part);
    }
}