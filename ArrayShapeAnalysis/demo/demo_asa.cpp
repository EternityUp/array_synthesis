
#include "array_shape_analysis.h"


int main(int argc, char *argv[]) {
    ArrayProps mic_array = {0};
#ifndef __SIMULATE_POSITIONS__
    printf("practical mic array positions\n");
    FILE *fp;
    const char *filename = argv[1];
    fp = fopen(filename, "r");

    if (NULL == fp) {
        printf("fopen error!\n");
        return -1;
    }

    mic_array.element_num = ArrayElements;
    ReadArrayInfoFromFile(mic_array.location, fp);
    GetCenteredArray(mic_array.location);
    fclose(fp);

#else
    printf("fabricated mic array positions\n");
    ReadArrayInfoFromFloatArray(CircularArray3DPos, mic_array.location);
    mic_array.is_linear =
        IsLinearArray(mic_array.location, mic_array.array_direction, mic_array.array_normal);
    if (!mic_array.is_linear)
        mic_array.is_planar =
        IsPlanarArray(mic_array.location, mic_array.array_normal);
    mic_array.is_tridimensional = !(mic_array.is_linear || mic_array.is_planar);

    ComputeElementDistancesCov(mic_array.location,mic_array.mic_spacings);
    mic_array.min_mic_spacing = GetMinimumMicSpacing(mic_array.mic_spacings);
#endif
    Point a = {1, 2, 3};
    Point b = {2, 3, 4};
    float c1 = PointNorm(a);
    float c2 = PointDistance(a, b);
    Point c = CrossProduct(a, b);
    float c3 = DotProduct(a, b);
    bool d1 = AreParallel(a, b);
    bool d2 = ArePerpendicular(a, b);
    bool d3 = AnyPoint(a);

    return 0;
}
