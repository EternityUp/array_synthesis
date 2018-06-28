#include "audio_record_process_thread.h"
#include "mic_array_processing.h"
#include "audio_file_process.h"
#include <pthread.h>
#include <sys/time.h>

extern int fd_cp;
extern int in_fd_proc;
extern int out_fd_proc;

const int kChannels = 6;

void* audio_capture(void *audio_capture_args)
{
  int ret;
  size_t write_size = 0, written = 0;
  audio_cp_data *capture_data = (audio_cp_data *)audio_capture_args;
  static size_t start_ind = 0;
  while(1)
  {
    
    ret = snd_pcm_readi(capture_data->pcm_handle, capture_data->audio_buffer, capture_data->num_frames);
    
    if (ret == -EPIPE)
    {
      fprintf(stderr, "overrun occured!\n");
      snd_pcm_prepare(capture_data->pcm_handle);
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
    else if (ret != (int)capture_data->num_frames)
    {
      fprintf(stderr, "short read, read %d frames\n", ret);
    }
    
    if (ret > 0)
    {
      memcpy(capture_data->cp_audio_buffer + start_ind,capture_data->audio_buffer,capture_data->audio_buffer_size);
      start_ind += capture_data->audio_buffer_size;
      if(start_ind == capture_data->cp_audio_buffer_size)
      {
	start_ind = 0;
	*capture_data->cond_state = 1;
      }
    }
    
    
    while(write_size != capture_data->audio_buffer_size)
    {
      written = write(fd_cp, capture_data->audio_buffer, capture_data->audio_buffer_size);
      if(written < 0)
      {
	fprintf(stderr,"write error in audio capture thread: %s\n",strerror(errno));
	fprintf(stderr,"exit from audio capture thread\n");
	pthread_exit(0);
      }
      write_size += written;
    }
    write_size = 0;
        
  }
  
}



void *audio_process(void *audio_process_args)
{
  size_t write_size = 0, written = 0;
  int write_size_out = 0, written_out = 0;
  audio_cp_data *proc_data = (audio_cp_data *)audio_process_args;
  
  // original input buffer for recorded audio stream    ------ 6 channels 
  const int original_multi_chs_frame_len = kChannels * kFrameLen;
  int16_t original_in_frame_buf_int16[original_multi_chs_frame_len] = {0};
  float original_in_frame_buf_float[original_multi_chs_frame_len] = {0.0f};
  float **original_in_multi_chs_buf = new float *[kChannels];
  for (int i = 0; i < kChannels; i++)
    original_in_multi_chs_buf[i] = new float[kFrameLen]{0.0};
  
  
  
  
  
  // input buffer for process   ------ 4 channels

  float **in_multi_chs_buf = new float *[kNumChannels];
  for (int i = 0; i < kNumChannels; i++)
    in_multi_chs_buf[i] = new float[kFrameLen]{0.0};
  
  // output buffer after process  ------ 1 channels
  float out_frame_buf_float[kFrameLen] = {0.0f};
  int16_t out_frame_buf_int16[kFrameLen] = {0};
  
  // variable for timing 
  struct timeval start_val, end_val;
  float elapsed_time;
  
  // instance of MicArrayProcInst 
  MicArrayProcInst pProc;
  // init mic array processing instance 
  InitMicArrayProc(&pProc);
  pProc.mic_array_info->azimuth = 0.0f;
  InitMicArrayFromFloatArray(XmosArrayPos, pProc.mic_array_info);
  
  // set parameters 
  int win_type = 1;
  SetParasMicArrayProc(&pProc, win_type);

  // validate parameters 
  ValidateParasMicArrayProc(&pProc);
  
  
  
  
  while(1)
  {
    if(*(proc_data->cond_state))
    {
      
      *(proc_data->cond_state) = 0;
      
      while(write_size != proc_data->audio_buffer_size)
      {
	written = write(in_fd_proc, proc_data->audio_buffer, proc_data->audio_buffer_size);
	if(written < 0)
	{
	  fprintf(stderr,"write error in audio process thread: %s\n",strerror(errno));
	  fprintf(stderr,"exit from audio process thread\n");
	  // free mic array processing instance 
          FreeMicArrayProc(&pProc);
	  for (int i = 0; i < kNumChannels; i++)
	    delete[] in_multi_chs_buf[i];
	  delete[] in_multi_chs_buf;
	  for (int i = 0; i < kChannels; i++)
	    delete[] original_in_multi_chs_buf[i];
	  delete[] original_in_multi_chs_buf;
	  pthread_exit(0);
	}
	write_size += written;
      }
      write_size = 0;
      
      // variable for file read and write operations 
      gettimeofday(&start_val, NULL);
      memcpy(original_in_frame_buf_int16,proc_data->audio_buffer,proc_data->audio_buffer_size);
      Int16_tToFloat(original_in_frame_buf_int16, original_multi_chs_frame_len, original_in_frame_buf_float);
      Deinterleave(original_in_frame_buf_float, original_in_multi_chs_buf, kChannels, kFrameLen);
      
      for (int i = 0; i < kNumChannels; i++)
      {
	memcpy(in_multi_chs_buf[i],original_in_multi_chs_buf[i+1],kFrameLen*sizeof(float));
      }
	

      // perform array processing algorithms 
      ProcCoreMicArrayProcSingleOut(&pProc, original_in_multi_chs_buf, out_frame_buf_float);
      // write processed audio data into output file 
      FloatToInt16_t(out_frame_buf_float, kFrameLen, out_frame_buf_int16);
      
      while(write_size_out != kFrameLen * sizeof(int16_t))
      {
	written_out = write(out_fd_proc, out_frame_buf_int16, kFrameLen * sizeof(int16_t));
	if(written_out < 0)
	{
	  fprintf(stderr,"write error in audio process thread: %s\n",strerror(errno));
	  fprintf(stderr,"exit from audio process thread\n");
	  // free mic array processing instance 
          FreeMicArrayProc(&pProc);
	  for (int i = 0; i < kNumChannels; i++)
	    delete[] in_multi_chs_buf[i];
	  delete[] in_multi_chs_buf;
	  for (int i = 0; i < kChannels; i++)
	    delete[] original_in_multi_chs_buf[i];
	  delete[] original_in_multi_chs_buf;
	  pthread_exit(0);
	}
	write_size_out += written_out;
      }
      write_size_out = 0;
      
      gettimeofday(&end_val, NULL);
      elapsed_time = (int64_t)((end_val.tv_sec - start_val.tv_sec) * 1000) + 1.f * (end_val.tv_usec - start_val.tv_usec) / 1000.f;
      
      fprintf(stderr,"Wasting time of each frame is :%.4f(ms)\n",elapsed_time);
      
     

    }
  }
}
