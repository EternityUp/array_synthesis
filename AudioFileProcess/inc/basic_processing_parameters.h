#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_BASIC_PROCESSING_PARAMETERS_H_
#define AMY_MODULES_AUDIO_PROCESSING_BASIC_PROCESSING_PARAMETERS_H_

#include "array_info.h"
#include "eigen_wrapper.h"



enum { kNumChannels = ArrayElements }; /* number of audio channels                       */
enum { kSampleRate = 16000 };          /* sample frequency                               */
enum { kFrameTime = 10 };              /* time duration per frame: ms                    */
enum { kFrameLen = 160 };              /* number of sample points per frame              */
enum { kBufferLen = 384 };             /* length of audio ring buffer to store frames    */
enum { kBlockLen = 256 };              /* length of audio block for tf analysis          */
enum { kBlockShift = 128 };            /* number of shift points of adjacent blocks      */
enum { kBlockLapped = 128 };           /* number of lapped points of adjacent blocks     */
enum { kBlockInitialDelay = 224 };     /* number of foregoing zero points of first block */
enum { kNumFreqs = 129 };              /* number of frequency bins                       */
enum { NIS = 10 };                     /* preceding frames used for bg estimation        */
enum { kMaxBlocks = 2 };               /* number of frequency bins                       */



typedef struct BlockTdFdInfo {
    MatrixXf in_td_block;     /* 2D: kBlockLen * kNumChannels */
    MatrixXf in_fd_block;     /* 2D: kNumFreqs * kNumChannels */
    MatrixXf ang_in_fd_block; /* 2D: kNumFreqs * kNumChannels */
    MatrixXf out_td_block;    /* 2D: kBlockLen * kNumChannels */
    MatrixXf out_fd_block;    /* 2D: kNumFreqs * kNumChannels */
} BlockTFInfo;


typedef struct DataBufferInfo {
    ColVectorXf in_buffer;       /* 1D : kBufferLen                                   */
    ColVectorXf out_buffer;      /* 1D : kBufferLen                                   */
    int read_pos = 160;          /* location of reader pointer in buffer              */
    int write_pos = 0;           /* lcoation of writer pointer in buffer              */
    bool rw_wrap = true;         /* relative location status of reader/writer pointer */
} BufferRWInfo;


#endif