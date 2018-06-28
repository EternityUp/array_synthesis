/*
 * record real-time audio via xmos usb sound card with multiple available channels,
 * and process them based on mic array processing algorithms.
 * audio-record and audio-process are divided into two threads.
 * so there are three threads:
 * one for main thread, namely manager thread
 * one for audio recording
 * and another one for processing 
 */

#include <alsa/asoundlib.h>
#include <iostream>
#include <fcntl.h>
#include "alsa_audio_record.h"
#include "audio_record_process_thread.h"
#include "basic_processing_parameters.h"


using std::cout;
using std::endl;

int fd_cp;   // file descriptor for capture audio stream
int in_fd_proc; // file descriptor for input process audio stream
int out_fd_proc; // file descriptor for output process audio stream

int main(int argc, char **argv) {
  #define __func__ "main function"
    fprintf(stderr, "----------------------------\n");
    fprintf(stderr, "start of alsa real-time capture and process demo!\n");
    
    if(argc < 2)
    {
      cout << "You must specify the sound card device!!!" << endl;
      return -1;
    }
    
    fd_cp = open("./record.pcm",O_WRONLY | O_CREAT | O_TRUNC,S_IRWXU);
    in_fd_proc = open("./in_process.pcm",O_WRONLY | O_CREAT | O_TRUNC,S_IRWXU);
    out_fd_proc = open("./out_process.pcm",O_WRONLY | O_CREAT | O_TRUNC,S_IRWXU);
    
    snd_pcm_t *capture_hdl = 0;
    snd_pcm_hw_params_t *capture_params;
    char *capture_buffer = 0, *proc_buffer = 0;
    ulong capture_buf_size, proc_buffer_size;
    
    char *capture_device = argv[1];
    uint sample_rate = 16000;
    uint num_channels = 6;
    snd_pcm_uframes_t capture_num_frames = 32;
    const uint ByteWidth = sizeof(int16_t);
    
    fprintf(stderr,"information about audio device and stream:\n \
    capture_device_name = %s \n \
    sample_rate = %d \n \
    num_channels = %d \n \
    capture_num_frames = %lu \n",
    capture_device,sample_rate,num_channels,capture_num_frames);
    
    pthread_t capture_thread, process_thread;
    int cond_state = 0;
    audio_cp_data cap_audio_cp_data, process_audio_cp_data;
    int err, *ret_cap_thread, *ret_proc_thread;
    
    // initialize capture device  
    alsa_device_init(capture_device,&capture_hdl,&capture_params,&sample_rate,&capture_num_frames,num_channels,SND_PCM_STREAM_CAPTURE);
    
    // allocate a buffer large enough to hold one period for capture device 
    alsa_alloc_buffer(&capture_buffer, &capture_buf_size, num_channels, capture_num_frames);
    memset(capture_buffer,0,capture_buf_size);
    
    // allocate a buffer for audio-process thread, its buffer size must be no smaller than that of capture buffer 
    proc_buffer_size = kFrameLen*num_channels*ByteWidth;
    proc_buffer = new char[proc_buffer_size];
    
    
    
    // create thread for real-time audio capture
    cap_audio_cp_data.pcm_handle = capture_hdl;
    cap_audio_cp_data.audio_buffer = capture_buffer;
    cap_audio_cp_data.cp_audio_buffer = proc_buffer;
    cap_audio_cp_data.cp_audio_buffer_size = proc_buffer_size;
    cap_audio_cp_data.num_frames = capture_num_frames;
    cap_audio_cp_data.audio_buffer_size = capture_buf_size;
    cap_audio_cp_data.cond_state = &cond_state;
    
    pthread_create(&capture_thread,NULL,audio_capture, (void*)(&cap_audio_cp_data));
    
    // create thread for real-time audio process
    process_audio_cp_data.pcm_handle = 0;
    process_audio_cp_data.audio_buffer = proc_buffer;
    process_audio_cp_data.cp_audio_buffer = 0;
    process_audio_cp_data.cp_audio_buffer_size = 0;
    process_audio_cp_data.num_frames = kFrameLen;
    process_audio_cp_data.audio_buffer_size = proc_buffer_size;
    process_audio_cp_data.cond_state = &cond_state;
       
    pthread_create(&process_thread,NULL,audio_process, (void*)(&process_audio_cp_data));
    
    // join audio-capture thread
    err = pthread_join(capture_thread, (void**)&ret_cap_thread);
    if(err)
      fprintf(stderr,"%s: %s\n", __func__, strerror(err));
    
    // join audio-process thread
    err = pthread_join(process_thread, (void**)&ret_proc_thread);
    if(err)
      fprintf(stderr,"%s: %s\n", __func__, strerror(err));
    
    
    
    // release resources
    alsa_release(capture_hdl,capture_buffer);
    delete[] proc_buffer;
    
    fprintf(stderr, "----------------------------\n");
    fprintf(stderr, "end of alsa real-time capture and process demo!\n");
    
  
    return 0;
}
