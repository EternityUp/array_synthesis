#pragma once

#include "alsa_audio_record.h"


typedef struct audio_capture_or_playback_data
{
  snd_pcm_t *pcm_handle;
  char *audio_buffer;
  char *cp_audio_buffer;
  snd_pcm_uframes_t num_frames;
  ulong audio_buffer_size;
  ulong cp_audio_buffer_size;
  int *cond_state;
}audio_cp_data;




void* audio_capture(void *audio_capture_args);


void *audio_process(void *audio_process_args);




