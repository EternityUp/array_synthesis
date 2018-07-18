#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_FFT_TRANSFORM_H
#define AMY_MODULES_AUDIO_PROCESSING_FFT_TRANSFORM_H

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

#include <complex>
#include "fft4g.h"

typedef std::complex<float> Complexf;


AMY_AUDIO_API bool IsPowerOf2(size_t n);

AMY_AUDIO_API size_t PowernOf2(size_t n);

AMY_AUDIO_API size_t NextPower2(size_t n);

AMY_AUDIO_API void FftTransform(const float *const in_td_data, Complexf *out_fd_data, size_t len);

AMY_AUDIO_API void InvFftTransform(const Complexf *const in_fd_data, float *out_td_data, size_t len);

AMY_AUDIO_API void ConstructComplexSpectrum(Complexf *complex_spectrum, float *power_spectrum,
                                            float *phase_spectrum, int nfreq_bins);

#endif