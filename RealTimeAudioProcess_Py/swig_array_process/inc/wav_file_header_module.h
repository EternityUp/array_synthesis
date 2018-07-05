#pragma once
#ifndef WAV_FILE_HEADER_MODULE_H_
#define WAV_FILE_HEADER_MODULE_H_

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


#include <stddef.h>
#include <stdint.h>
#include <limits>


typedef std::numeric_limits<int16_t> limits_int16;
static const size_t kWavHeaderSize = 44;

class AMY_AUDIO_API ReadableWav {
public:
    // Returns the number of bytes read.
    size_t virtual Read(void *buf, size_t num_bytes) = 0;

    virtual ~ReadableWav() {}
};

enum WavFormat {
    kWavFormatPcm = 1,  // PCM, each sample of size bytes_per_sample
    kWavFormatALaw = 6,  // 8-bit ITU-T G.711 A-law
    kWavFormatMuLaw = 7,  // 8-bit ITU-T G.711 mu-law
};

// Return true if the given parameters will make a welled WAV header.
bool CheckWavParameters(size_t num_channels,
                        int samle_rate,
                        WavFormat format,
                        size_t bytes_per_sample,
                        size_t num_samples);

// Write a kWavHeaderSize bytes long WAV header to buf. The payload that 
// follows the header is supposed to have the specified number of interleaved
// channels and contain the specified total number of samples of the specified
// type. Checks the input parameters for validity.
void WriteWavHeader(uint8_t *buf,
                    size_t num_channels,
                    int sample_rate,
                    WavFormat format,
                    size_t bytes_per_sample,
                    size_t num_samples);

// Read a WAV header from an implemented ReadableWav and parse the values into
// the provided output parameters. ReadableWav is used because the header can
// be variably sized. Returns false if the header is invalid.
bool ReadWavHeader(ReadableWav *readable,
                   size_t *num_channels,
                   int *sample_rate,
                   WavFormat *format,
                   size_t *bytes_per_sample,
                   size_t *num_samples);


void FloatS16ToS16(const float *src, size_t size, int16_t *dest);

static inline int16_t FloatS16ToS16(float v) {
    static const float kMaxRound = std::numeric_limits<int16_t>::max() - 0.5f;
    static const float kMinRound = std::numeric_limits<int16_t>::min() + 0.5f;
    if (v > 0)
        return v >= kMaxRound ? std::numeric_limits<int16_t>::max()
                              : static_cast<int16_t>(v + 0.5f);
    return v <= kMinRound ? std::numeric_limits<int16_t>::min() : static_cast<int16_t>(v - 0.5f);
}


#endif
