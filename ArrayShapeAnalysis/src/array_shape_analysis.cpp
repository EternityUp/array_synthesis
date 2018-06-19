#include "array_shape_analysis.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

const float kMaxDotProduct = 1e-6f;

AMY_AUDIO_API void GetCenteredArray(Point *mic_array) {
    Point centered_point = {0};
    for (size_t i = 0; i < ArrayElements; i++) {
        centered_point.x += mic_array[i].x;
        centered_point.y += mic_array[i].y;
        centered_point.z += mic_array[i].z;
    }
    centered_point.x /= ArrayElements;
    centered_point.y /= ArrayElements;
    centered_point.z /= ArrayElements;
    for (size_t i = 0; i < ArrayElements; i++) {
        mic_array[i].x -= centered_point.x;
        mic_array[i].y -= centered_point.y;
        mic_array[i].z -= centered_point.z;
    }
}


AMY_AUDIO_API bool AnyPoint(const Point &a) {
    if (fabs(a.x) > kMaxDotProduct || fabs(a.y) > kMaxDotProduct || fabs(a.z) > kMaxDotProduct)
        return true;
    return false;
}

AMY_AUDIO_API Point PairDirection(const Point &a, const Point &b) {
    return {b.x - a.x, b.y - a.y, b.z - a.z};
}

AMY_AUDIO_API float DotProduct(const Point &a, const Point &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

AMY_AUDIO_API Point CrossProduct(const Point &a, const Point &b) {
    Point tmp;
    tmp.x = a.y * b.z - a.z * b.y;
    tmp.y = a.z * b.x - a.x * b.z;
    tmp.z = a.x * b.y - a.y * b.x;
    return tmp;
}


AMY_AUDIO_API bool AreParallel(const Point &a, const Point &b) {
    Point tmp = CrossProduct(a, b);
    return AnyPoint(tmp) ? false : true;
}

AMY_AUDIO_API bool ArePerpendicular(const Point &a, const Point &b) {
    return fabs(DotProduct(a, b)) > kMaxDotProduct ? false : true;
}

AMY_AUDIO_API int IsLinearArray(const Point *const mic_array, Point &direction, Point &normal) {
    Point first_pair_direction = PairDirection(mic_array[0], mic_array[1]);
    for (size_t i = 2; i < ArrayElements; i++) {
        Point pair_direction = PairDirection(mic_array[i - 1], mic_array[i]);
        if (!AreParallel(first_pair_direction, pair_direction))
            return 0;
    }
    direction = first_pair_direction;
    normal = {-first_pair_direction.y, first_pair_direction.x, 0};
    return 1;
}

AMY_AUDIO_API int IsPlanarArray(const Point *const mic_array, Point &normal) {
    if (ArrayElements < 3)
        return 0;

    Point first_pair_direction = PairDirection(mic_array[0], mic_array[1]);
    Point second_pair_direction = PairDirection(mic_array[1], mic_array[2]);
    Point tmp_normal = CrossProduct(first_pair_direction, second_pair_direction);

    for (size_t i = 3; i < ArrayElements; i++) {
        Point other_pair_direction = PairDirection(mic_array[i - 1], mic_array[i]);
        if (!ArePerpendicular(tmp_normal, other_pair_direction))
            return 0;
    }
    normal = tmp_normal;
    return 1;
}


AMY_AUDIO_API float PointNorm(const Point &a) {
    return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}

AMY_AUDIO_API float PointDistance(const Point &a, const Point &b) {
    Point tmp = {a.x - b.x, a.y - b.y, a.z - b.z};
    return PointNorm(tmp);
}

AMY_AUDIO_API void ComputeElementDistancesCov(const Point *const mic_array, float c[][ArrayElements]) {
    for (size_t i = 0; i < ArrayElements; i++)
        for (size_t j = 0; j < ArrayElements; j++) {
            *(*(c + i) + j) = PointDistance(mic_array[i], mic_array[j]);
        }
}

AMY_AUDIO_API float GetMinimumMicSpacing(float(*c)[ArrayElements]) {
    float min_spacing = 100.f;
    for (size_t i = 0; i < ArrayElements; i++)
        for (size_t j = i + 1; j < ArrayElements; j++) {
            float cij = *(*(c + i) + j);
            min_spacing = min_spacing < cij ? min_spacing : cij;
        }
    return min_spacing;
}

AMY_AUDIO_API float ComputeSoundPath(const Point &a, float angle_rad) {
    return a.x * cosf(angle_rad) + a.y * sinf(angle_rad);
}

AMY_AUDIO_API void ComputeElemetsSoundPath(const Point *const mic_array, float angle_rad, float *sound_path) {
    for (int i = 0; i < ArrayElements; i++) {
        sound_path[i] = ComputeSoundPath(mic_array[i], angle_rad);
    }
}

AMY_AUDIO_API void ReadArrayInfoFromFile(Point *mic_array, FILE *fp) {
    const int MAX_LINE = 1024;
    const int BITS = 4;
    const int HEADLINES = 3;
    char buf[MAX_LINE];
    int len = 0, line_num = 0;
    while (NULL != fgets(buf, MAX_LINE, fp)) {
        len = strlen(buf);
        buf[len - 1] = '\0';
        line_num++;
        printf("%d %s\n", line_num, buf);
        if (line_num > HEADLINES && line_num <= HEADLINES + ArrayElements) {
            char *decimal_point, buf_float[BITS] = {'\0'};
            int i = 2;
            while (NULL != (decimal_point = strrchr(buf, '.'))) {
                memcpy(buf_float, decimal_point - 1, BITS);
                switch (i) {
                    case 2:
                        mic_array[line_num - HEADLINES - 1].z = float(atof(buf_float));
                        break;
                    case 1:
                        mic_array[line_num - HEADLINES - 1].y = float(atof(buf_float));
                        break;
                    case 0:
                        mic_array[line_num - HEADLINES - 1].x = float(atof(buf_float));
                        break;
                    default:
                        break;

                }
                i--;
                *decimal_point = '\0';
            }
        }
    }
}

AMY_AUDIO_API void ReadArrayInfoFromFloatArray(const float array1[ArrayElements][3], Point *array2) {
    for (size_t i = 0; i < ArrayElements; i++) {
        array2[i].x = array1[i][0];
        array2[i].y = array1[i][1];
        array2[i].z = array1[i][2];
    }
}


AMY_AUDIO_API void InitMicArrayFromTxtFile(const char *filename, ArrayProps *pMicArray) {

    FILE *array_fp;

    array_fp = fopen(filename, "r");

    if (array_fp == NULL) {
        printf("Failure of opening array information file: %s\n ", filename);
    }

    pMicArray->element_num = ArrayElements;
    ReadArrayInfoFromFile(pMicArray->location, array_fp);
    GetCenteredArray(pMicArray->location);
    pMicArray->is_linear =
            IsLinearArray(pMicArray->location, pMicArray->array_direction, pMicArray->array_normal);
    if (!pMicArray->is_linear)
        pMicArray->is_planar =
                IsPlanarArray(pMicArray->location, pMicArray->array_normal);
    pMicArray->is_tridimensional = !(pMicArray->is_linear || pMicArray->is_planar);
    ComputeElementDistancesCov(pMicArray->location, pMicArray->mic_spacings);
    pMicArray->min_mic_spacing = GetMinimumMicSpacing(pMicArray->mic_spacings);
    ComputeElemetsSoundPath(pMicArray->location, pMicArray->azimuth, pMicArray->sound_path);

    fclose(array_fp);

}

AMY_AUDIO_API void InitMicArrayFromFloatArray(const float array1[ArrayElements][3], ArrayProps *pMicArray) {
    pMicArray->element_num = ArrayElements;
    ReadArrayInfoFromFloatArray(array1, pMicArray->location);
    GetCenteredArray(pMicArray->location);
    pMicArray->is_linear =
            IsLinearArray(pMicArray->location, pMicArray->array_direction, pMicArray->array_normal);
    if (!pMicArray->is_linear)
        pMicArray->is_planar =
                IsPlanarArray(pMicArray->location, pMicArray->array_normal);
    pMicArray->is_tridimensional = !(pMicArray->is_linear || pMicArray->is_planar);
    ComputeElementDistancesCov(pMicArray->location, pMicArray->mic_spacings);
    pMicArray->min_mic_spacing = GetMinimumMicSpacing(pMicArray->mic_spacings);
    ComputeElemetsSoundPath(pMicArray->location, pMicArray->azimuth, pMicArray->sound_path);

}