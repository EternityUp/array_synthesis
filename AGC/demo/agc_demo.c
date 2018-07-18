#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "agc_wrapper.h"


/*
int main(int argc, char *argv[])
{
	FILE *in_fp = NULL;
	FILE *out_fp = NULL;
	void *agcHandle = NULL;	
	short *pData    = NULL;
	short *pOutData = NULL;

	int fs = 16000;
	int minLevel = 0;
	int maxLevel = 255;
	int agcMode  = kAgcModeAdaptiveAnalog;
	WebRtcAgc_config_t agcConfig;  

	int frameSize = 160;
	int len = frameSize*sizeof(short);
	int micLevelIn = 0;
	int micLevelOut = 0;
	uint8_t saturationWarning;
	int nAgcRet = 0;
	int read_bytes = 0;
	int write_bytes = 0;
	int frame_num = 0;


	in_fp =  fopen(argv[1], "rb");
	out_fp = fopen(argv[2], "wb");

	if(NULL == in_fp || NULL == out_fp)
	{
		printf("fopen error!\n");
		system("pause");
		return -1;
	}


	

	WebRtcAgc_Create(&agcHandle);
	WebRtcAgc_Init(agcHandle, minLevel, maxLevel, agcMode, fs);

	
	agcConfig.compressionGaindB = 15;
	agcConfig.limiterEnable     = 1;
	agcConfig.targetLevelDbfs   = 3;
	WebRtcAgc_set_config(agcHandle, agcConfig);


	pData    = (short*)malloc(frameSize*sizeof(short));
	pOutData = (short*)malloc(frameSize*sizeof(short));  

	
	while(!feof(in_fp))
	{
		memset(pData, 0, len);
		memset(pOutData, 0, len);
		read_bytes = fread(pData, 1, len, in_fp);
		if(read_bytes < 0)
		{
			printf("fread error!\n");
			system("pause");
			return -1;
		}
		//nAgcRet = WebRtcAgc_AddMic(agcHandle, pData, NULL, frameSize);
		//nAgcRet = WebRtcAgc_VirtualMic(agcHandle, pData, NULL, frameSize, micLevelIn, &micLevelOut);
		nAgcRet = WebRtcAgc_Process(agcHandle, pData, NULL, frameSize, pOutData,NULL, micLevelIn, &micLevelOut, 0, &saturationWarning);
		if (nAgcRet == -1)
		{
			printf("WebRtcAgc_Process error!\n");
			system("pause");
			return -1;

		}
		micLevelIn = micLevelOut;
		write_bytes = fwrite(pOutData, 1, len, out_fp);
		if (write_bytes < 0)
		{
			printf("fwrite error!\n");
			system("pause");
			return -1;
		}
		frame_num++;
		printf("i = %d, micLevelOut=%d\n", frame_num, micLevelOut);
	}

	fclose(in_fp);
	fclose(out_fp);
	free(pData);
	free(pOutData);

	WebRtcAgc_Free(agcHandle);

	system("pause");

	return 0;
}
*/


int main(int argc, char *argv[])
{
	FILE *in_fp = NULL;
	FILE *out_fp = NULL;
	short *pData = NULL;
	short *pOutData = NULL;
	Agc_Outer_Config agc_config;

	agc_config.agc_inst = NULL;
	agc_config.sample_frequency = 16000;
	agc_config.min_mic_level = 0;
	agc_config.max_mic_level = 255;
	agc_config.in_mic_level = 0;
	agc_config.out_mic_level = 0;
	agc_config.agc_mode = kAgcModeAdaptiveDigital;
	agc_config.compressionGaindB = 20;
	agc_config.targetLevelDbfs = 6;
	agc_config.limiterEnable = 1;
	agc_config.frame_len = 160;

	int len = agc_config.frame_len * sizeof(short);
	int nAgcRet = 0;
	int read_bytes = 0;
	int write_bytes = 0;
	int frame_num = 0;

	in_fp = fopen(argv[1], "rb");
	out_fp = fopen(argv[2], "wb");

	if (NULL == in_fp || NULL == out_fp)
	{
		printf("fopen error!\n");
		system("pause");
		return -1;
	}

	agc_Create(&agc_config);
	agc_Init(&agc_config);
	agc_Configure(&agc_config);

	pData = (short*)malloc(agc_config.frame_len * sizeof(short));
	pOutData = (short*)malloc(agc_config.frame_len * sizeof(short));


	while (!feof(in_fp))
	{
		memset(pData, 0, len);
		memset(pOutData, 0, len);
		read_bytes = fread(pData, 1, len, in_fp);
		if (read_bytes < 0)
		{
			printf("fread error!\n");
			system("pause");
			return -1;
		}
		nAgcRet = agc_ProcessCore(pData, pOutData, &agc_config);
		if (nAgcRet == -1)
		{
			printf("WebRtcAgc_Process error!\n");
			system("pause");
			return -1;

		}
		write_bytes = fwrite(pOutData, 1, len, out_fp);
		if (write_bytes < 0)
		{
			printf("fwrite error!\n");
			system("pause");
			return -1;
		}
		frame_num++;
	}

	fclose(in_fp);
	fclose(out_fp);
	free(pData);
	free(pOutData);

	agc_Free(&agc_config);


	return 0;
}