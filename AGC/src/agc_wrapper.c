#include "agc_wrapper.h"

AMY_AUDIO_API int agc_Create(Agc_Outer_Config * agc_config)
{
	return WebRtcAgc_Create(&(agc_config->agc_inst));
}

AMY_AUDIO_API int agc_Init(Agc_Outer_Config * agc_config)
{
	int32_t minLevel, maxLevel;
	uint32_t fs;
	int16_t agc_mode;
	minLevel = agc_config->min_mic_level;
	maxLevel = agc_config->max_mic_level;
	agc_mode = agc_config->agc_mode;
	fs = agc_config->sample_frequency;
	return WebRtcAgc_Init(agc_config->agc_inst, minLevel, maxLevel, agc_mode, fs);
}

AMY_AUDIO_API int agc_Configure(Agc_Outer_Config * agc_config)
{
	WebRtcAgc_config_t agcConfig;
	agcConfig.compressionGaindB = agc_config->compressionGaindB;
	agcConfig.targetLevelDbfs = agc_config->targetLevelDbfs;
	agcConfig.limiterEnable = agc_config->limiterEnable;
	return WebRtcAgc_set_config(agc_config->agc_inst, agcConfig);
}

AMY_AUDIO_API int agc_ProcessCore(const int16_t *in_data,
	int16_t *out_data, Agc_Outer_Config * agc_config)
{
	uint8_t saturationWarning = 0;

	return WebRtcAgc_Process(agc_config->agc_inst, in_data, 0, agc_config->frame_len, out_data,
		0, agc_config->in_mic_level, &agc_config->out_mic_level, 0, &saturationWarning);
}

AMY_AUDIO_API int agc_Free(Agc_Outer_Config * agc_config)
{
	return WebRtcAgc_Free(agc_config->agc_inst);
}

AMY_AUDIO_API void agc_CvtToInt16(float *src_floatS16, int16_t *dst, int len)
{
	for (int i = 0; i < len; i++)
	{
		dst[i] = (int16_t)(src_floatS16[i] + 0.5f);
	}
}

AMY_AUDIO_API void agc_CvtInt16ToFloatS16(int16_t *src, float *dst_floatS16, int len)
{
	for (int i = 0; i < len; i++)
	{
		dst_floatS16[i] = src[i] * 1.0f;
	}
}
