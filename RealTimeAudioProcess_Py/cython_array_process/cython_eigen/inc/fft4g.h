#pragma once
#ifndef COMMON_AUDIO_FFT4G_H_
#define COMMON_AUDIO_FFT4G_H_

#if defined(__cplusplus)
extern "C" {
#endif

// Refer to fft4g.c for documentation.
void rdft(size_t n, int isgn, float *a, size_t *ip, float *w);

#if defined(__cplusplus)
}
#endif

#endif  // COMMON_AUDIO_FFT4G_H_
