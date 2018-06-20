#include "basic_signal_processing.h"

#define LEN 6

int main(int argc, char *argv[]) {
    int length = 16;
    float *Hanning_win = new float[length];
    float *Kaiser_win = new float[length];
    float alpha = 1.5f;

    Hanning(length, Hanning_win);
    KaiserBesselDerived(alpha, length, Kaiser_win);

    delete[] Hanning_win;
    delete[] Kaiser_win;

    int16_t vector1[LEN] = {1, 2, 3, 2, 1, 4};
    int16_t vector2[LEN] = {0, 2, 1, 1, 1, 3};
    int scaling = 1;
    int32_t dpws = DotProductWithScale(vector1, vector2, LEN, scaling);

    float src[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    float dst[10] = {0};
    float win[3] = {0};
    Hanning(3, win);
    Smooth1Ddata(src, dst, 10, 3, win);

    float norm = VectorNorm(src, 10);
    ScaleVector(src, dst, 10, 1.f / norm);

    size_t gcd = ComputeGcd(100, 24);

    
    return 0;
}