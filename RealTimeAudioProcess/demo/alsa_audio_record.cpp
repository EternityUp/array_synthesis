#include "alsa_audio_record.h"
#include <stdio.h>


void alsa_device_init(const char *device_name,
		       snd_pcm_t **pcm_handle,
		       snd_pcm_hw_params_t **params,
		       uint *sample_rate,
		       snd_pcm_uframes_t *num_frames,
		       uint num_channels,
		       snd_pcm_stream_t stream)
{
  int dir = 0;
  int ret = 0;
  
  // open pcm device 
  int mode = 0;
  ret = snd_pcm_open(pcm_handle,device_name,stream,mode);
  if(ret < 0)
  {
    fprintf(stderr, "Unable to open pcm device:%s\n err = %s\n",\
    device_name, snd_strerror(ret));
    exit(OPEN_ERROR);
  }
 
  // allocate a hardware parameters object 
  snd_pcm_hw_params_alloca(params);
  
  // fill hardware parameters with default value
  snd_pcm_hw_params_any(*pcm_handle, *params);
  
  // set desired hardware parameters
  // Interleaved mode
  snd_pcm_hw_params_set_access(*pcm_handle,*params,SND_PCM_ACCESS_RW_INTERLEAVED);
  
  // Signed 16-bit little-endian format
  snd_pcm_hw_params_set_format(*pcm_handle,*params,SND_PCM_FORMAT_S16_LE);
  
  // channels 
  snd_pcm_hw_params_set_channels(*pcm_handle,*params,num_channels);
  
  // sample rate 
  snd_pcm_hw_params_set_rate_near(*pcm_handle,*params,(unsigned int *)sample_rate,&dir);
  
  // period size
  snd_pcm_hw_params_set_period_size_near(*pcm_handle,*params,(snd_pcm_uframes_t*)num_frames,&dir);
  snd_pcm_hw_params_get_period_size(*params, (snd_pcm_uframes_t*)num_frames, &dir);
  
  // Write the parameters to the driver
  ret = snd_pcm_hw_params(*pcm_handle,*params);
  if(ret < 0)
  {
    fprintf(stderr,
            "unable to set hw parameters: %s\n",
            snd_strerror(ret));
    exit(SET_ERROR);
  }
  
}




void alsa_alloc_buffer(char **buffer,ulong *buf_size,uint num_channels,snd_pcm_uframes_t num_frames)
{

  *buf_size = num_frames * num_channels * 2; 
  *buffer = (char *)malloc(*buf_size);
}



void alsa_print_params(snd_pcm_hw_params_t *params)
{
  ;
}



void alsa_capture_process(snd_pcm_t *capture_hdl, char *capture_buf, snd_pcm_uframes_t capture_frames, ulong capture_buf_size)
{
  int ret = 0;
  ret = snd_pcm_readi(capture_hdl, capture_buf, capture_frames);
  
  if (ret == -EPIPE)
  {
    fprintf(stderr, "overrun occured!\n");
    snd_pcm_prepare(capture_hdl);
  }
  else if (ret == -EBADFD)
  {
    fprintf(stderr, "pcm device is not in the right state!\n");
  }
  else if (ret == -ESTRPIPE)
  {
    fprintf(stderr,"stream is suspended and waiting for an application recovery!\n");
  }
  else if (ret < 0)
  {
    fprintf(stderr, "snd_pcm_readi error: %s\n", snd_strerror(ret));
  }
  else if (ret != (int)capture_frames)
  {
    fprintf(stderr, "short read, read %d frames\n", ret);
  }
  
}



void alsa_release(snd_pcm_t *pcm_handle, char *buffer)
{
  snd_pcm_drain(pcm_handle);
  snd_pcm_close(pcm_handle);
  free(buffer);
}

