#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_ARRAY_INFO_H_
#define AMY_MODULES_AUDIO_PROCESSING_ARRAY_INFO_H_

#define M_PI_F  3.14159265358979323846f   // pi

const int ArrayElements = 4;


//const float XmosArrayPos[ArrayElements][3] = {
//	{ 0.0430f, 0.0750f, 0.0000f },   /* 实际对应板上mic1*/
//{ 0.0430f, 0.0000f, 0.0000f },       /* 实际对应板上mic3*/
//{ 0.0000f, 0.0000f, 0.0000f },       /* 实际对应板上mic4*/
//{ 0.0000f, 0.0750f, 0.0000f },       /* 实际对应板上mic6*/
//};


/* 自行测量值 */
const float XmosArrayPos[ArrayElements][3] = {
	{ 0.04337f, 0.0000f, 0.00000f },
{ -0.02168f, -0.03756f, 0.00000f },
{ -0.04337f, 0.00000f, 0.00000f },
{ 0.02168f, 0.03756f, 0.00000f },
};



//const float CircularArray3DPos[ArrayElements][3] = {
//	{ 2.2000f, 3.2000f, 0.5000f },
//{ 2.1500f, 3.2866f, 0.5000f },
//{ 2.0500f, 3.2866f, 0.5000f },
//{ 2.0000f, 3.2000f, 0.5000f },
//{ 3.0500f, 3.1134f, 0.5000f },
//{ 2.1500f, 3.1134f, 0.5000f }
//};
//
//
//const float LinearArray3DPos[ArrayElements][3] = {
//	{ 0.5000f, 2.5000f, 3.0000f },
//{ 0.5500f, 2.6000f, 3.1500f },
//{ 0.6000f, 2.7000f, 3.3000f },
//{ 0.6500f, 2.8000f, 3.4500f },
//{ 0.7000f, 2.9000f, 3.6000f },
//{ 0.7500f, 3.0000f, 3.7500f }
//};
//
//const float ArbitraryArray3D[ArrayElements][3] = {
//	{ 0.5377f,     -0.4336f,    0.7254f },
//{ 1.8339f,     0.3426f,     -0.0631f },
//{ -2.2588f,    3.5784f,     0.7147f },
//{ 0.8622f,     2.7694f,     -0.2050f },
//{ 0.3188f,     -1.3499f,    -0.1241f },
//{ -1.3077f,    3.0349f,     1.4897f }
//};

const float radian_resolution = M_PI_F / 180;



#endif
