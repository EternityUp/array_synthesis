#pragma once 


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
#include <stdint.h>
#include "gain_control.h"
	typedef struct agc_outer_config
	{
		uint32_t sample_frequency;
		int16_t agc_mode;
		int32_t min_mic_level;
		int32_t max_mic_level;
		int32_t in_mic_level;
		int32_t out_mic_level;
		int16_t targetLevelDbfs;   // default 3 (-3 dBOv)
		int16_t compressionGaindB; // default 9 dB
		uint8_t limiterEnable;     // default kAgcTrue (on)
		int16_t frame_len;
		void *agc_inst;
	}Agc_Outer_Config;


	AMY_AUDIO_API  int agc_Create(Agc_Outer_Config * agc_config);

	AMY_AUDIO_API  int agc_Init(Agc_Outer_Config * agc_config);

	AMY_AUDIO_API  int agc_Configure(Agc_Outer_Config * agc_config);

	AMY_AUDIO_API  int agc_ProcessCore(const int16_t *in_data, 
		int16_t *out_data, Agc_Outer_Config * agc_config);

	AMY_AUDIO_API  int agc_Free(Agc_Outer_Config * agc_config);

	AMY_AUDIO_API void agc_CvtToInt16(float *src_floatS16, int16_t *dst, int len);

	AMY_AUDIO_API void agc_CvtInt16ToFloatS16(int16_t *src, float *dst_floatS16, int len);


#ifdef __cplusplus
}
#endif /* __cplusplus*/