#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_ARRAY_SHAPE_ANALYSIS_H_
#define AMY_MODULES_AUDIO_PROCESSING_ARRAY_SHAPE_ANALYSIS_H_

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

#include <stdio.h>

#include "array_info.h"

#define MAX_LINE_LENGTH 1024

/* Cartesian coordinate system */
typedef struct CartsianPoint {
    float x;
    float y;
    float z;
} Point;

/* Spherical coordinate system */
typedef struct SphericalPoint {
    float azimuth;
    float elevation;
    float radius;
} SphericalPointf;


/* 3D properties of mic array */
typedef struct ArrayProperties3D {
    Point location[ArrayElements];                     /* array element coordinates in 3D cartesian space        */
    Point array_normal;                                /* normal direction of an array if it is linear or planar */
    Point array_direction;                             /* direction of an array if it is linear                  */
    float min_mic_spacing;                             /* minimum spacing between array elements                 */
    float mic_spacings[ArrayElements][ArrayElements];  /* spacings between array elements                        */
    bool is_linear;
    bool is_planar;
    bool is_tridimensional;
    int element_num;
    float azimuth;
    float sound_path[ArrayElements];                   /* sound path of every elements for centered mic array    */
} ArrayProps;


AMY_AUDIO_API void GetCenteredArray(Point *mic_array);

AMY_AUDIO_API bool AnyPoint(const Point &a);

AMY_AUDIO_API Point PairDirection(const Point &a, const Point &b);

AMY_AUDIO_API float DotProduct(const Point &a, const Point &b);

AMY_AUDIO_API Point CrossProduct(const Point &a, const Point &b);

AMY_AUDIO_API bool AreParallel(const Point &a, const Point &b);

AMY_AUDIO_API bool ArePerpendicular(const Point &a, const Point &b);

AMY_AUDIO_API int IsLinearArray(const Point *const mic_array, Point &direction, Point &normal);

AMY_AUDIO_API int IsPlanarArray(const Point *const mic_array, Point &normal);

AMY_AUDIO_API float PointNorm(const Point &a);

AMY_AUDIO_API float PointDistance(const Point &a, const Point &b);

AMY_AUDIO_API void ComputeElementDistancesCov(const Point *const mic_array, float c[][ArrayElements]);

AMY_AUDIO_API float GetMinimumMicSpacing(float(*c)[ArrayElements]);

AMY_AUDIO_API float ComputeSoundPath(const Point &a, float angle_rad);

AMY_AUDIO_API void ComputeElemetsSoundPath(const Point *const mic_array, float angle_rad, float *sound_path);

AMY_AUDIO_API inline float DegreesToRadians(float angle_degrees) {
    return M_PI_F * angle_degrees;
}

AMY_AUDIO_API inline float RadiansToDegrees(float angle_radians) {
    return 180.0f * angle_radians / M_PI_F;
}

AMY_AUDIO_API void ReadArrayInfoFromFile(Point *mic_array, FILE *fp);

AMY_AUDIO_API void ReadArrayInfoFromFloatArray(const float array1[ArrayElements][3], Point *array2);

AMY_AUDIO_API void InitMicArrayFromTxtFile(const char *filename, ArrayProps *pMicArray);

AMY_AUDIO_API void InitMicArrayFromFloatArray(const float array1[ArrayElements][3], ArrayProps *pMicArray);


#ifdef __cplusplus
}
#endif /* __cplusplus*/

#endif