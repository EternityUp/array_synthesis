#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_MIC_ARRAY_PROCESSING_H_
#define AMY_MODULES_AUDIO_PROCESSING_MIC_ARRAY_PROCESSING_H_

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

#include "basic_processing_parameters.h"
#include "array_shape_analysis.h"
#include "fft_transform.h"


typedef struct VoiceAndNoiseEstInst {
    const int16_t smoothed_win_len_of_sp = 3;       /* size of smoothed window for spectrum */
    const float alpha_sp = 0.8f;                    /* smoothing factor for spectrum  */
    const int16_t noise_minsp_est_win_len = 125;    /* window size for minsp est      */
    const float likelihood_ratio_threshold = 1.0f;  /* threshold for speech decision  */
    const float alpha_speech_presence_prob = 0.8f;  /* smoothing factor for speech presence probability */
    const float beta_prior_snr = 0.7f;              /* smoothing factor for prior snr estimation */
    const int16_t smoothed_local_win_len = 3;       /* size of local smoothing window for average prior snr */
    const int16_t smoothed_global_win_len = 31;     /* size of global smoothing window for average prior snr */
    const float compromising = 0.99f;               /* weighting factor of noise suppression and speech distortion*/
    const float sigma_min = 0.316228f;
    const float sigma_max = 0.562341f;
    const float sigma_pmin = 1.000000f;
    const float sigma_pmax = 3.162278f;
    const float max_speech_absence_prob = 0.95f;    /* maximum prior speech absence probability */
    const float Gmin = 0.0316228f;                  /* minimum attenuation factor for noise */
    const float prior_speech_presence = 0.01f;
    const float kKbdAlpha = 1.5f;
    const float alpha_ns_prob_update = 0.95f;       /* smoothing factor for updating noise spectrum estimation */
    const float snr_prev_val = 0.0562341f;

    ColVectorXf *smooth_win_of_sp;                   /* smoothing window for spectrum                 */
    ColVectorXf *smoothed_local_win;                 /* local smoothing window for average prior snr  */
    ColVectorXf *smoothed_global_win;                /* global smoothing window for average prior snr */

} VaNEstInst;


typedef struct VoiceActivityDetectionInst {

    /* processed frequency range */
    const int16_t frequency_start = 50;
    const int16_t frequency_end = 4000;
    int16_t frequency_start_ind = 0;
    int16_t frequency_end_ind = 0;


    /* voice decision threshold */
    const float beyond_ratio_energy = 1.0f;       /* energy statis                  */
    const float beyond_ratio_ee = 1.0f;           /* energy-to-entropy ratio statis */
    const float initial_voice_prob = 0.05f;       /* initial voice probability      */
    const float voice_prob_threshold = 0.2f;      /* voice probability              */

    /* smoothing factor */
    const float smoothing = 0.80f;
    const int16_t smoothed_win_len_of_sp = 3;     /* size of smoothed window for statistics */
    ColVectorXf *smooth_win_of_sp;                /* smoothing window for statistics        */
    /* number of leading larger statistic value */
    const int LargerNum = 10;

    const float alpha_pre_emphasis = 0.9375f;     /* pre_emphasis factor                 */

    const float alpha_log_sp = 2.0f;              /* factor used to compute log spectrum */


    /* speech segment decision */
    float speech_start_threshold = 0.2f;
    float speech_finished_threshold = 0.2f;
    const uint16_t non_speech_to_speech_frame_len = 20;
    const uint16_t speech_to_non_speech_frame_len = 20;

    /* initialize speech property to be false */
    bool is_speech_frame = false;
    bool is_speech_triggered = false;

} VoiceDetInst;


typedef struct DirectionOfArrivalEstInst {
    /* processed frequency range */
    const float frequency_start = 50;
    const float frequency_end = 4000;
    int16_t frequency_start_ind = 0;
    int16_t frequency_end_ind = 0;
    float freq_bin[kNumFreqs] = {0.0f};

    /* information of rough sweep */
    float coarse_ang_span = 10.0f;           /* angle span for rough sweep                           */
    float coarse_ang_start = 0.0f;           /* initial angle for rough sweep                        */
    float coarse_ang_max1 = 180.0f;          /* terminated angle for rough sewwp (line array)        */
    float coarse_ang_max2 = 360.0f;          /* terminated angle for rough sewwp (planar array)      */
    float coarse_beam_num = kNumChannels;    /* to be considered (rough + fine sweep)                */

    /* information of fine sweep */
    float fine_center_ang = 0.0f;            /* centered angle of fine sweep                         */
    float fine_ang_span = 3.0f;              /* angle span for fine sweep                            */

    /* loading factor for adaptive beamforming */
    float loading_factor = 0.01f;

    float sp_smoothing = 0.8f;               /* smoothing factor for complex spectrum                */

    /* doa estimation result */
    float doa = 0.0f;


} DoaInst;


typedef struct ArrayBeamformingInst {
    float freq_bin[kNumFreqs] = {0.0f};
} bfInst;


typedef struct PostFilterInst {
    const int Num_Inters = kNumChannels * (kNumChannels - 1) / 2;
    const float smooth_alpha = 0.8f;
    const float min_postfilter_coeff = 0.05f;
    const float freq_low = 200.0f;
    const float freq_high = 6800.0f;
    int freq_low_ind = 0;
    int freq_high_ind = 0;
    const float slope = 0.7f;
    ColVectorXf *min_postfilter_coeff_cn;
} PfInst;


typedef struct MIC_ARRAY_PROC_INST {
    /* Structural pointer for information of mic array */
    ArrayProps *mic_array_info;

    /* Structural pointer for estimation of speech and noise, further, snr */
    VaNEstInst *snr_est;

    /* Structural pointer for voice frame determination */
    VoiceDetInst *vadInst;
    int vad_mode = 0;    /*
						 0: speech probability 
						 1: energy
						 2: energy-entropy ratio 
						 */


    /* Structural pointer for source direction estimation */
    DoaInst *doa_core;
    int doa_mode = 0;   /* 0:cbf; 1:mvdr; 2:stmv */


    /* Structural pointer for beamforming */
    bfInst *bf_core;
    int bf_mode = 0; /* 0:cbf */

    /* Structural pointer for postfilter */
    PfInst *pf_core;
    int cn_mode = 0;


    /* Structural variable for input and output audio ring buffer */
    BufferRWInfo buf_rw_status[kNumChannels];

    /* Structural variable for decomposed block data of every frame data */
    BlockTFInfo block_tf[kMaxBlocks];

    /* window to process time domain block data */
    ColVectorXf window;

    /* output of speech and noise spectrum estimation module */
    MatrixXf noise_sps_est;
    MatrixXf recur_avg_sps_est;
    MatrixXf prior_snr;
    MatrixXf post_snr;
    MatrixXf gain;
    MatrixXf speech_prob;

    /* voice detection result */
    bool is_speech = false;
    bool is_speech_triggered = false;

    /* doa estimation result */
    float doa = 0.0f;


    uint32_t frame_num = 0;                       /* sequence number of current frame      */
    int frame_offset = 0;                         /* offset of audio ring buffer           */
    const float kKbdAlpha = 1.5f;                 /* parameter for computing Kaiser window */
    float freq_resolution = 0.0f;                 /* frequency resolution                  */

} MicArrayProcInst;


AMY_AUDIO_API void InitMicArrayProc(MicArrayProcInst *self);

AMY_AUDIO_API void SetParasMicArrayProc(MicArrayProcInst *self, int win_type);

AMY_AUDIO_API void ValidateParasMicArrayProc(MicArrayProcInst *self);

AMY_AUDIO_API void
ProcCoreMicArrayProcSingleOut(MicArrayProcInst *self, float **in_multi_chs_data, float *out_single_ch_data);

AMY_AUDIO_API void
ProcCoreMicArrayProcMultiOut(MicArrayProcInst *self, float **in_multi_chs_data, float **out_multi_chs_data);

AMY_AUDIO_API void ComputePowerDensitySp(float *sp, const Complexf *xf, int len);

AMY_AUDIO_API void ComputeAngularSp(float *ang, const Complexf *xf, int len);

AMY_AUDIO_API void
EstimateSnrAndGainEtc(MatrixXf *xf, VaNEstInst *psnr, MatrixXf *noise_est, MatrixXf *recur_avg_sp_est,
                      MatrixXf *prior_snr, MatrixXf *post_snr, MatrixXf *gain, MatrixXf *speech_prob, int frame_seq,
                      int ch_seq);

AMY_AUDIO_API void RecursiveAverage(float *prev, float *curr, float *dst, int len, float factor);

AMY_AUDIO_API void TranverseMinOf2Vecs(float *vec1, float *vec2, float *dst, int len);

AMY_AUDIO_API void TranverseMaxOf2Vecs(float *vec1, float *vec2, float *dst, int len);

AMY_AUDIO_API void ComputeDotQuotientOf2Vecs(float *num_vec, float *den_vec, float *dst, int len);

AMY_AUDIO_API void ComputeDotProductOf2Vecs(float *vec1, float *vec2, float *dst, int len);

AMY_AUDIO_API void CompareVecWithScalar(float *vec1, float scalar, float *dst, int len);

AMY_AUDIO_API void RecurAveVecWithScalar(float *vec, float scalar, float *dst, int len, float factor);

AMY_AUDIO_API void RecursiveAverageWithSmoothVec(float *prev, float *curr, float *dst, int len, float *factor);

AMY_AUDIO_API void CopyValsOfVecsToOther(float *src, float *dst, int len);

AMY_AUDIO_API void ConstrainLocalGlobalSnr(float *local_snr, float *global_snr, float min_th, float max_th, int len);

AMY_AUDIO_API float ComputeVecHeadAvg(float *vec, int start_ind, int end_ind);

AMY_AUDIO_API void ComputePFrameAndSnrPeak(float *p_frame, float *snr_peak, float *snr_peak_prev, float seg_avg_snr,
                                           float *seg_avg_snr_prev, float snr_min, float snr_max, float snr_pmin,
                                           float snr_pmax, int ch_seq);

AMY_AUDIO_API void ComputeSpeechAbsenceProb(float *absence_prob, float *local_snr,
                                            float *global_snr, float frame_snr, int len);

AMY_AUDIO_API void GetMinOfVecAndScalar(float *vec, float *min_vec, float scalar, int len);

AMY_AUDIO_API void GetMaxOfVecAndScalar(float *vec, float *max_vec, float scalar, int len);

AMY_AUDIO_API void ComputePriorSNRCurr(float *prior_snr, float *gain_prev,
                                       float *post_snr_prev, float *post_snr_curr, float compromising, int len);

AMY_AUDIO_API void ComputeVkPost_Prior(double *vk, float *prior_snr_curr, float *post_snr_curr, int len);

AMY_AUDIO_API void ComputeGain_Curr(float *gain_curr, float *prior_snr_curr, double *exp_integeral, int len);

AMY_AUDIO_API void ComputeSpeechPresenceProb(float *speech_absence_prob, float *prior_snr_curr,
                                             double *vk_post_prior, float *speech_presence_prob, int len);

AMY_AUDIO_API void ComputeFinalGain(float *final_gain, float *gain_curr,
                                    float *speech_presence_prob, float gain_min, int len);

AMY_AUDIO_API void ApplyGainToSp(float *gain, float *original_sp, float *processed_sp, int len);

AMY_AUDIO_API void VoiceDecisionWithinFrame(VoiceDetInst *vad, float *avg_sf, float *avg_nf,
                                            float *voice_prob, int frame_seq, int vad_mode);

AMY_AUDIO_API void ValidVoiceSegDecisionAlongFrames(VoiceDetInst *vad, int frame_seq);

AMY_AUDIO_API int ComputeSumOfBoolVec(int *bvec, int len);

AMY_AUDIO_API void DoaEstimationCore(MatrixXcf *chs_com_sp, DoaInst *pdoa,
                                     ArrayProps *mic_info, int doa_mode);

AMY_AUDIO_API void SetComplexfVec2Zeros(Complexf(*vec2D)[kNumChannels], int rows);

AMY_AUDIO_API void DoaEstComputeSteerVec(float *sound_path, Complexf *steer_vcetor, float frequency);

AMY_AUDIO_API Complexf ComputeCompVecConjDotProduct(const Complexf *vec1, const Complexf *vec2, int len);

AMY_AUDIO_API float
ComputeOutSpatialSpectrumCbf(const Complexf *steer_vcetor, const Complexf *In_FreqSpectrum, int len);

AMY_AUDIO_API float NormComplex(Complexf &a);

AMY_AUDIO_API Complexf ComputeCompVecDotProduct(const Complexf *vec1, const Complexf *vec2, int len);

AMY_AUDIO_API float
ComputeOutSpatialSpectrumMvdr(const Complexf *steer_vcetor, const Complexf *In_FreqSpectrum, int len,
                              float load_factor);

AMY_AUDIO_API void ComputeStmvNormCov(const Complexf *vec1, const Complexf *vec2, Complexf(*norm_mat)[kNumChannels]);


AMY_AUDIO_API float ComputeOutSpatialSpectrumStmv(Complexf(*norm_mat)[kNumChannels], int len, float load_factor);


AMY_AUDIO_API void BfComputeSteerVec(float *sound_path, Complexf *steer_vcetor, float frequency);

AMY_AUDIO_API void
BfCore(MatrixXcf *in_chs_com_sp, MatrixXcf *out_chs_com_sp, bfInst *pbf, ArrayProps *mic_info, int bf_mode);

AMY_AUDIO_API void BfPhaseCompensation(const Complexf *steer_vec, Complexf *src_sp, Complexf *dst_sp, int len);

AMY_AUDIO_API void PostFilterCore(MatrixXcf *in_chs_com_sp, MatrixXcf *single_out_ch_com_sp, PfInst *ppf, int cn_mode);

AMY_AUDIO_API void
PostFilterGetFilterCoeffMin(float *min_filter_coeff, int low_ind, int high_ind, int len, float slope);

AMY_AUDIO_API void ComputePowerDensity(float *power_density_prev, float *power_density,
                                       Complexf *ch1_comp_spectrum, Complexf *ch2_comp_spectrum, float smooth_factor,
                                       int len);

AMY_AUDIO_API void ComputeAveragePowerDensity(float *src_power_density, float *dst_power_density, int chs, int freqs);

AMY_AUDIO_API void ComputePostFilterCoeff(float *voice_sp, float *noisy_voice_sp, float *filter_coeff, int len);

AMY_AUDIO_API void GetCwiseMaxValOfVecAndScalar(float *vec1, float scalar, float *dst, int len);

AMY_AUDIO_API void GetCwiseMinValOf2Vecs(float *vec1, float *vec2, float *dst, int len);

AMY_AUDIO_API void GetCwiseMaxValOf2Vecs(float *vec1, float *vec2, float *dst, int len);

AMY_AUDIO_API void GetCwiseMinValOfVecAndScalar(float *vec1, float scalar, float *dst, int len);

AMY_AUDIO_API void
ComputeAverageComSpectrum(Complexf *src_com_spectrum, Complexf *dst_com_spectrum, int chs, int freqs);

AMY_AUDIO_API void ApplyPostFilterCoeff(Complexf *src_sp, Complexf *dst_sp, float *filter_coeff, int len);

AMY_AUDIO_API void ConstructComSpFromEsp_Psp(Complexf *dst_sp, float *energy_sp, float *phase_sp, int len);

AMY_AUDIO_API void FreeMicArrayProc(MicArrayProcInst *self);


#endif


