#include "basic_signal_processing.h"
#include <math.h>
#include <complex>
#include <assert.h>
#include <string.h>

using std::complex;
static const long double epsilon = 10.0 * LDBL_EPSILON;

// Modified Bessel function of order 0 for complex inputs.
static complex<float> I0(complex<float> x) {
	complex<float> y = x / 3.75f;
	y *= y;
	return 1.0f + y * (
		3.5156229f + y * (
			3.0899424f + y * (
				1.2067492f + y * (
					0.2659732f + y * (
						0.360768e-1f + y * 0.45813e-2f)))));
}


AMY_AUDIO_API void Hanning(int length, float* window) {
	assert(length >= 1);
	assert(window != nullptr);
	for (int i = 1; i <= length; ++i) {
		window[i - 1] = 0.5f * (1 - cosf(2 * static_cast<float>(M_PI) * i /
			(length + 1)));
	}
}


AMY_AUDIO_API void KaiserBesselDerived(float alpha, size_t length,
	float* window) {
	assert(length >= 1);
	assert(window != nullptr);

	const size_t half = (length + 1) / 2;
	float sum = 0.0f;

	for (size_t i = 0; i <= half; ++i) {
		complex<float> r = (4.0f * i) / length - 1.0f;
		sum += I0(static_cast<float>(M_PI) * alpha * sqrt(1.0f - r * r)).real();
		window[i] = sum;
	}
	for (size_t i = length - 1; i >= half; --i) {
		window[length - i - 1] = sqrtf(window[length - i - 1] / sum);
		window[i] = window[length - i - 1];
	}
	if (length % 2 == 1) {
		window[half - 1] = sqrtf(window[half - 1] / sum);
	}
}


AMY_AUDIO_API int32_t DotProductWithScale(const int16_t* vector1,
	const int16_t* vector2,
	size_t length,
	int scaling)
{
	int64_t sum = 0;
	size_t i = 0;

	/* Unroll the loop to improve performance. */
	for (i = 0; i + 3 < length; i += 4) {
		sum += (vector1[i + 0] * vector2[i + 0]) >> scaling;
		sum += (vector1[i + 1] * vector2[i + 1]) >> scaling;
		sum += (vector1[i + 2] * vector2[i + 2]) >> scaling;
		sum += (vector1[i + 3] * vector2[i + 3]) >> scaling;
	}
	for (; i < length; i++) {
		sum += (vector1[i] * vector2[i]) >> scaling;
	}

	return int32_t(sum);
}

AMY_AUDIO_API void Smooth1Ddata(float *src, float *dst, int data_len, int win_len, const float *window)
{
	assert(win_len % 2 == 1 && win_len > 1);
	int half_win_len = (win_len - 1) / 2;
	int start_ind = 0, end_ind = 0;
	for (int i = 0; i < data_len; i++)
	{
		dst[i] = 0; // must be initialized to be zero
		start_ind = 0 > i - half_win_len ? 0 : i - half_win_len;
		end_ind = data_len - 1 < i + half_win_len ? data_len - 1 : i + half_win_len;
		for (int j = start_ind; j <= end_ind; j++)
			dst[i] += src[j] * window[j - start_ind];
		dst[i] /= (end_ind - start_ind + 1);
	}
}


AMY_AUDIO_API void ScaleVector(float *src, float *dst, float length, float scale)
{
	for (int i = 0; i < length; i++)
	{
		dst[i] = src[i] * scale;
	}
}

AMY_AUDIO_API float VectorNorm(float *src, float length)
{
	float tmp = 0.0f;
	for (int i = 0; i < length; i++)
	{
		tmp += src[i] * src[i];
	}
	return sqrtf(tmp);
}

AMY_AUDIO_API size_t ComputeGcd(size_t a, size_t b)
{
	assert(a > 0 && b > 0);
	size_t temp;
	while (b)
	{
		temp = a;
		a = b;
		b = temp % b;
	}
	return a;
}


AMY_AUDIO_API void ApplyWindow1D(float *src, float *dst, float *window, int win_len)
{
	for (int i = 0; i < win_len; i++)
	{
		dst[i] = src[i] * window[i];
	}
}

AMY_AUDIO_API void ApplyWindow2D(float **src, int rows, int cols, float *window)
{
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			src[i][j] *= window[j];
		}
}


AMY_AUDIO_API double Exponential_Integral_Ei(double x)
{
	return (double)xExponential_Integral_Ei((long double)x);
}

AMY_AUDIO_API long double xExponential_Integral_Ei(long double x)
{
	if (x < -5.0L) return Continued_Fraction_Ei(x);
	if (x == 0.0L) return -DBL_MAX;
	if (x < 6.8L)  return Power_Series_Ei(x);
	if (x < 50.0L) return Argument_Addition_Series_Ei(x);
	return Continued_Fraction_Ei(x);
}


AMY_AUDIO_API long double Continued_Fraction_Ei(long double x)
 {
	 long double Am1 = 1.0L;
	 long double A0 = 0.0L;
	 long double Bm1 = 0.0L;
	 long double B0 = 1.0L;
	 long double a = expl(x);
	 long double b = -x + 1.0L;
	 long double Ap1 = b * A0 + a * Am1;
	 long double Bp1 = b * B0 + a * Bm1;
	 int j = 1;

	 a = 1.0L;
	 while (fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1)) {
		 if (fabsl(Bp1) > 1.0L) {
			 Am1 = A0 / Bp1;
			 A0 = Ap1 / Bp1;
			 Bm1 = B0 / Bp1;
			 B0 = 1.0L;
		 }
		 else {
			 Am1 = A0;
			 A0 = Ap1;
			 Bm1 = B0;
			 B0 = Bp1;
		 }
		 a = -j * j;
		 b += 2.0L;
		 Ap1 = b * A0 + a * Am1;
		 Bp1 = b * B0 + a * Bm1;
		 j += 1;
	 }
	 return (-Ap1 / Bp1);
 }


AMY_AUDIO_API long double Power_Series_Ei(long double x)
 {
	 long double xn = -x;
	 long double Sn = -x;
	 long double Sm1 = 0.0L;
	 long double hsum = 1.0L;
	 long double g = 0.5772156649015328606065121L;
	 long double y = 1.0L;
	 long double factorial = 1.0L;

	 if (x == 0.0L) return (long double)-DBL_MAX;

	 while (fabsl(Sn - Sm1) > epsilon * fabsl(Sm1)) {
		 Sm1 = Sn;
		 y += 1.0L;
		 xn *= (-x);
		 factorial *= y;
		 hsum += (1.0 / y);
		 Sn += hsum * xn / factorial;
	 }
	 return (g + logl(fabsl(x)) - expl(x) * Sn);
 }

AMY_AUDIO_API long double Argument_Addition_Series_Ei(long double x)
 {
	 static long double ei[] = {
		 1.915047433355013959531e2L,  4.403798995348382689974e2L,
		 1.037878290717089587658e3L,  2.492228976241877759138e3L,
		 6.071406374098611507965e3L,  1.495953266639752885229e4L,
		 3.719768849068903560439e4L,  9.319251363396537129882e4L,
		 2.349558524907683035782e5L,  5.955609986708370018502e5L,
		 1.516637894042516884433e6L,  3.877904330597443502996e6L,
		 9.950907251046844760026e6L,  2.561565266405658882048e7L,
		 6.612718635548492136250e7L,  1.711446713003636684975e8L,
		 4.439663698302712208698e8L,  1.154115391849182948287e9L,
		 3.005950906525548689841e9L,  7.842940991898186370453e9L,
		 2.049649711988081236484e10L, 5.364511859231469415605e10L,
		 1.405991957584069047340e11L, 3.689732094072741970640e11L,
		 9.694555759683939661662e11L, 2.550043566357786926147e12L,
		 6.714640184076497558707e12L, 1.769803724411626854310e13L,
		 4.669055014466159544500e13L, 1.232852079912097685431e14L,
		 3.257988998672263996790e14L, 8.616388199965786544948e14L,
		 2.280446200301902595341e15L, 6.039718263611241578359e15L,
		 1.600664914324504111070e16L, 4.244796092136850759368e16L,
		 1.126348290166966760275e17L, 2.990444718632336675058e17L,
		 7.943916035704453771510e17L, 2.111342388647824195000e18L,
		 5.614329680810343111535e18L, 1.493630213112993142255e19L,
		 3.975442747903744836007e19L, 1.058563689713169096306e20L
	 };
	 int  k = (int)(x + 0.5);
	 int  j = 0;
	 long double xx = (long double)k;
	 long double dx = x - xx;
	 long double xxj = xx;
	 long double edx = expl(dx);
	 long double Sm = 1.0L;
	 long double Sn = (edx - 1.0L) / xxj;
	 long double term = DBL_MAX;
	 long double factorial = 1.0L;
	 long double dxj = 1.0L;

	 while (fabsl(term) > epsilon * fabsl(Sn)) {
		 j++;
		 factorial *= (long double)j;
		 xxj *= xx;
		 dxj *= (-dx);
		 Sm += (dxj / factorial);
		 term = (factorial * (edx * Sm - 1.0L)) / xxj;
		 Sn += term;
	 }

	 return ei[k - 7] + Sn * expl(xx);
 }

AMY_AUDIO_API void Exponential_Integral_Ei_vector(double *x, double *y, int len)
{
	for (int i = 0; i < len; i++)
	{
		y[i] = Exponential_Integral_Ei(x[i]);
	}
}



AMY_AUDIO_API void PreEmphasis(float *src, float *dst, float alpha, int len)
{
	float *tmp = new float[len] {0};
	dst[0] = src[0];
	tmp[0] = src[0];
	for (int i = 1; i < len; i++)
	{
		tmp[i] = src[i];
		dst[i] = tmp[i] - alpha * tmp[i - 1];
	}

	delete[] tmp;
}


AMY_AUDIO_API float ComputeIntervalAverage(float *src, int range_ind_start, int range_ind_end)
{
	float tmp = 0.0f;
	for (int i = range_ind_start; i <= range_ind_end; i++)
		tmp += src[i];
	tmp = tmp / (range_ind_end - range_ind_start + 1);
	return tmp;
}


AMY_AUDIO_API float ComputeLogEnergyEntropyRatio(float *src, int range_ind_start, int range_ind_end, float alpha)
{
	float sum_energy = 0.0f, entropy = 0.0f, log_spectrum = 0.0f;
	int len = range_ind_end - range_ind_start + 1, i = 0;
	float *energy_pdf = new float[len] {0};
	for (i = range_ind_start; i <= range_ind_end; i++)
	{
		sum_energy += src[i];
	}
	int ret = range_ind_start;
	for (i = 0; i < len; i++)
	{
		energy_pdf[i] = src[ret] / sum_energy;
		if (energy_pdf[i] < 1e-6f)
			energy_pdf[i] = 1e-6f;
		ret++;
	}
	/* compute entropy */
	for (i = 0; i < len; i++)
	{
		entropy += -energy_pdf[i] * log2f(energy_pdf[i]);
	}

	delete[] energy_pdf;

	/* compute log spectrum */
	log_spectrum = log10f(1 + sum_energy / alpha);

	/* compute log energy to entropy ratio */
	return log_spectrum / entropy;

}

AMY_AUDIO_API void RecursiveSmoothing(float *src1, float*src2, float *dst, float alpha, int len)
{
	for (int i = 0; i < len; i++)
		dst[i] = alpha * src1[i] + (1 - alpha)*src2[i];
}