
#pragma once 

#include <alsa/asoundlib.h>


typedef unsigned int uint;
typedef unsigned long ulong;

typedef enum ERROR_CODE {
  OPEN_ERROR = 1,   /* open sound card device error */
  SET_ERROR  = 2    /* set device parameters error  */
}Error_Code;

#ifdef __cplusplus
extern "C" {
#endif

void alsa_device_init(const char *device_name,
		       snd_pcm_t **pcm_handle,
		       snd_pcm_hw_params_t **params,
		       uint *sample_rate,
		       snd_pcm_uframes_t *num_frames,
		       uint num_channels,
		       snd_pcm_stream_t stream
 		    );

void alsa_alloc_buffer(char **buffer,ulong *buf_size,uint num_channels,snd_pcm_uframes_t num_frames);

void alsa_print_params(snd_pcm_hw_params_t *params);

void alsa_capture_process(snd_pcm_t *capture_hdl, char *capture_buf, snd_pcm_uframes_t capture_frames, ulong capture_buf_size);

void alsa_release(snd_pcm_t *pcm_handle, char *buffer);

#ifdef __cplusplus
}
#endif