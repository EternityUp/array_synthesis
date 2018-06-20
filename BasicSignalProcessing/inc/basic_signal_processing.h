#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_BASIC_SIGNAL_PROCESSING_H_
#define AMY_MODULES_AUDIO_PROCESSING_BASIC_SIGNAL_PROCESSING_H_

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

#include <stdlib.h>

AMY_AUDIO_API void Hanning(int length, float *window);

AMY_AUDIO_API void KaiserBesselDerived(float alpha, size_t length, float *window);

AMY_AUDIO_API int32_t DotProductWithScale(const int16_t *vector1,
                                          const int16_t *vector2,
                                          size_t length,
                                          int scaling);

AMY_AUDIO_API void Smooth1Ddata(float *src, float *dst, int data_len, int win_len, const float *window);

AMY_AUDIO_API void ScaleVector(float *src, float *dst, float length, float scale);

AMY_AUDIO_API float VectorNorm(float *src, float length);

AMY_AUDIO_API size_t ComputeGcd(size_t a, size_t b);

AMY_AUDIO_API void ApplyWindow1D(float *src, float *dst, float *window, int win_len);

AMY_AUDIO_API void ApplyWindow2D(float **src, int rows, int cols, float *window);

AMY_AUDIO_API double Exponential_Integral_Ei(double x);

AMY_AUDIO_API long double xExponential_Integral_Ei(long double x);

AMY_AUDIO_API long double Continued_Fraction_Ei(long double x);

AMY_AUDIO_API long double Power_Series_Ei(long double x);

AMY_AUDIO_API long double Argument_Addition_Series_Ei(long double x);

AMY_AUDIO_API void Exponential_Integral_Ei_vector(double *x, double *y, int len);

AMY_AUDIO_API void PreEmphasis(float *src, float *dst, float alpha, int len);

AMY_AUDIO_API float ComputeIntervalAverage(float *src, int range_ind_start, int range_ind_end);

AMY_AUDIO_API float ComputeLogEnergyEntropyRatio(float *src, int range_ind_start, int range_ind_end, float alpha);

AMY_AUDIO_API void RecursiveSmoothing(float *src1, float *src2, float *dst, float alpha, int len);


#ifdef __cplusplus
}
#endif /* __cplusplus*/


#endif






