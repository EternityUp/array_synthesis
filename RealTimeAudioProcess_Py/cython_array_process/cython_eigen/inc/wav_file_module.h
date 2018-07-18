#pragma once
#ifndef WAVE_FILE_MODULE_H
#define WAVE_FILE_MODULE_H


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

#include <string>

// Interface to provide access to WAV file parameters.
class AMY_AUDIO_API WavFile {
public:
    virtual ~WavFile() {}

    virtual int sample_rate() const = 0;

    virtual size_t num_channels() const = 0;

    virtual size_t num_samples() const = 0;

    // Return a human-readable string containing the audio format.
    std::string FormatAsString() const;
};

// Simple C++ class for writing 16-bit PCM WAV files. 
class AMY_AUDIO_API WavWriter final : public WavFile {
public:
    // Open a new WAV file for writing.
    WavWriter(const std::string &filename, int sample_rate, size_t num_channels);

    // Close the WAV file, after writing its header.
    ~WavWriter() override;

    // Write additional samples to the file. Each sample is in the range
    // [-32768, 32767], and there must be the previously specified number of
    // interleaved channels.
    void WriteSamples(const float *sample, size_t num_samples);

    void WriteSamples(const int16_t *samples, size_t num_samples);

    int sample_rate() const override;

    size_t num_channels() const override;

    size_t num_samples() const override;

private:
    void Close();

    const int sample_rate_;
    const size_t num_channels_;
    size_t num_samples_;
    FILE *file_handle_;
};


// Follows the conventions of WavWriter.
class AMY_AUDIO_API WavReader final : public WavFile {
public:
    // Opens an existing WAV file for reading.
    explicit WavReader(const std::string &filename);

    // Close the WAV file.
    ~WavReader() override;

    // Returns the number of samples read. If this is less than requested,
    // verifies that the end of the file was reached.
    size_t ReadSamples(size_t num_samples, float *samples);

    size_t ReadSamples(size_t num_samples, int16_t *samples);

    int sample_rate() const override;

    size_t num_channels() const override;

    size_t num_samples() const override;

private:
    void Close();

    int sample_rate_;
    size_t num_channels_;
    size_t num_samples_; // Total number of samples int the file.
    size_t num_samples_remaining_;
    FILE *file_handle_; // Input file, owned by this class.
};


#endif

