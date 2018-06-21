#include "mic_array_processing.h"
#include "audio_file_process.h"
#include "basic_signal_processing.h"


const float sound_speed = 340;                /* sound speed in air: m/s              */


/* comparison function*/
static int cmp(const void *a, const void *b) {
    float ca = *(float *) a;
    float cb = *(float *) b;
    return (cb > ca) ? 1 : -1;  /* descending order */
    //	return (ca > cb) ? 1 : -1;  /* ascending order */
}

AMY_AUDIO_API void InitMicArrayProc(MicArrayProcInst *self) {
    int i = 0;

    self->mic_array_info = new ArrayProps[sizeof(ArrayProps)];

    self->snr_est = new VaNEstInst[sizeof(VaNEstInst)];
    self->snr_est->smooth_win_of_sp = new ColVectorXf(self->snr_est->smoothed_win_len_of_sp, 1);
    self->snr_est->smoothed_local_win = new ColVectorXf(self->snr_est->smoothed_local_win_len, 1);
    self->snr_est->smoothed_global_win = new ColVectorXf(self->snr_est->smoothed_global_win_len, 1);
    self->snr_est->smooth_win_of_sp->setZero(self->snr_est->smoothed_win_len_of_sp, 1);
    self->snr_est->smoothed_local_win->setZero(self->snr_est->smoothed_local_win_len, 1);
    self->snr_est->smoothed_global_win->setZero(self->snr_est->smoothed_global_win_len, 1);

    self->vadInst = new VoiceDetInst[sizeof(VoiceDetInst)];
    self->vadInst->smooth_win_of_sp = new ColVectorXf(self->vadInst->smoothed_win_len_of_sp, 1);
    self->vadInst->smooth_win_of_sp->setZero(self->vadInst->smoothed_win_len_of_sp, 1);

    self->doa_core = new DoaInst[sizeof(DoaInst)];

    self->bf_core = new bfInst[sizeof(bfInst)];

    self->pf_core = new PfInst[sizeof(PfInst)];
    self->pf_core->min_postfilter_coeff_cn = new ColVectorXf(int(kNumFreqs), 1);
    self->pf_core->min_postfilter_coeff_cn->setZero(kNumFreqs, 1);


    for (i = 0; i < kNumChannels; i++) {
        self->buf_rw_status[i].in_buffer.setZero(kBufferLen, 1);
        self->buf_rw_status[i].out_buffer.setZero(kBufferLen, 1);
    }

    for (i = 0; i < kMaxBlocks; i++) {
        self->block_tf[i].in_td_block.setZero(kBlockLen, kNumChannels);
        self->block_tf[i].in_fd_block.setZero(kNumFreqs, kNumChannels);
        self->block_tf[i].ang_in_fd_block.setZero(kNumFreqs, kNumChannels);
        self->block_tf[i].out_td_block.setZero(kBlockLen, kNumChannels);
        self->block_tf[i].out_fd_block.setZero(kNumFreqs, kNumChannels);
    }

    self->window.setZero(kBlockLen, 1);

    self->noise_sps_est.setZero(kNumFreqs, kNumChannels);
    self->recur_avg_sps_est.setZero(kNumFreqs, kNumChannels);
    self->prior_snr.setZero(kNumFreqs, kNumChannels);
    self->post_snr.setZero(kNumFreqs, kNumChannels);
    self->gain.setZero(kNumFreqs, kNumChannels);
    self->speech_prob.setZero(kNumFreqs, kNumChannels);
}


AMY_AUDIO_API void SetParasMicArrayProc(MicArrayProcInst *self, int win_type) {
    float norm = 0.0f;

    if (win_type == 0) {
        Hanning(kBlockLen, self->window.data());
    } else if (win_type == 1) {
        KaiserBesselDerived(self->kKbdAlpha, kBlockLen, self->window.data());
    }

    Hanning(self->snr_est->smoothed_win_len_of_sp, self->snr_est->smooth_win_of_sp->data());
    norm = VectorNorm(self->snr_est->smooth_win_of_sp->data(), self->snr_est->smoothed_win_len_of_sp);
    ScaleVector(self->snr_est->smooth_win_of_sp->data(), self->snr_est->smooth_win_of_sp->data(),
                self->snr_est->smoothed_win_len_of_sp, 1.0f / norm);

    Hanning(self->snr_est->smoothed_local_win_len, self->snr_est->smoothed_local_win->data());
    norm = VectorNorm(self->snr_est->smoothed_local_win->data(), self->snr_est->smoothed_local_win_len);
    ScaleVector(self->snr_est->smoothed_local_win->data(), self->snr_est->smoothed_local_win->data(),
                self->snr_est->smoothed_local_win_len, 1.0f / norm);

    Hanning(self->snr_est->smoothed_global_win_len, self->snr_est->smoothed_global_win->data());
    norm = VectorNorm(self->snr_est->smoothed_global_win->data(), self->snr_est->smoothed_global_win_len);
    ScaleVector(self->snr_est->smoothed_global_win->data(), self->snr_est->smoothed_global_win->data(),
                self->snr_est->smoothed_global_win_len, 1.0f / norm);

    self->freq_resolution = 1.0f * kSampleRate / kBlockLen;
    self->vadInst->frequency_start_ind = int(ceil(self->vadInst->frequency_start / self->freq_resolution)) - 1;
    self->vadInst->frequency_end_ind = int(ceil(self->vadInst->frequency_end / self->freq_resolution)) - 1;
    Hanning(self->vadInst->smoothed_win_len_of_sp, self->vadInst->smooth_win_of_sp->data());
    norm = VectorNorm(self->vadInst->smooth_win_of_sp->data(), self->vadInst->smoothed_win_len_of_sp);
    ScaleVector(self->vadInst->smooth_win_of_sp->data(), self->vadInst->smooth_win_of_sp->data(),
                self->vadInst->smoothed_win_len_of_sp, 1.0f / norm);

    for (int i = 0; i < kNumFreqs; i++) {
        self->doa_core->freq_bin[i] = i * self->freq_resolution;
        self->bf_core->freq_bin[i] = i * self->freq_resolution;
    }

    self->doa_core->frequency_start_ind = int(ceil(self->doa_core->frequency_start / self->freq_resolution)) - 1;
    self->doa_core->frequency_end_ind = int(ceil(self->doa_core->frequency_end / self->freq_resolution)) - 1;

    self->pf_core->freq_low_ind = int(std::ceil(self->pf_core->freq_low / self->freq_resolution)) - 1;
    self->pf_core->freq_high_ind = int(std::ceil(self->pf_core->freq_high / self->freq_resolution)) - 1;
    if (self->cn_mode == 1) {
        PostFilterGetFilterCoeffMin(self->pf_core->min_postfilter_coeff_cn->data(), self->pf_core->freq_low_ind,
                                    self->pf_core->freq_high_ind, kNumFreqs, self->pf_core->slope);
    }
}


AMY_AUDIO_API void ValidateParasMicArrayProc(MicArrayProcInst *self) {
    assert(kSampleRate == 16000);
}


AMY_AUDIO_API void
ProcCoreMicArrayProcSingleOut(MicArrayProcInst *self, float **in_multi_chs_data, float *out_single_ch_data) {
    int i = 0, kBlocks = 0, first_frame_in_block = 0;

    self->is_speech = false;

    /* read input audio data into audio ring buffer */
    for (i = 0; i < kNumChannels; i++) {
        InputBufferWrite(in_multi_chs_data[i], self->buf_rw_status[i].in_buffer.data(), kFrameLen,
                         self->buf_rw_status[i].read_pos, self->buf_rw_status[i].write_pos,
                         self->buf_rw_status[i].rw_wrap, kBufferLen);
    }

    /*
    based on input audio ring buffer, first_frame_in_block, etc,
    divide each frame data into blocks with kBlockLapped samples in adjacent blocks
    for each channel, process each frame block by block
    output data will be delayed by some samples
    */

    /* temporary variable to store complex spectrum */
    MatrixXcf tmp_fd, tmp_fd_chs, tmp_fd_chs2, single_ch_out_sp, avg_energy_sp;
    static MatrixXcf tmp_fd_chs_prev = MatrixXcf::Zero(kNumFreqs, kNumChannels);
    tmp_fd.setZero(kNumFreqs, 1);
    tmp_fd_chs.setZero(kNumFreqs, kNumChannels);
    tmp_fd_chs2.setZero(kNumFreqs, kNumChannels);
    single_ch_out_sp.setZero(kNumFreqs, 1);
    avg_energy_sp.setZero(kNumFreqs, 1);

    MatrixXf post_filter_output_block = MatrixXf::Zero(kBlockLen, 1);
    static MatrixXf post_filter_output_buffer = MatrixXf::Zero(kBufferLen, 1);


    /* temporary variable to store average pds used for voice detection */
    ColVectorXf chs_avg_sf, chs_avg_nf, chs_avg_speech_prob;

    first_frame_in_block = self->frame_offset;
    while (first_frame_in_block < kFrameLen) {
        chs_avg_sf.setZero(kNumFreqs, 1);
        chs_avg_nf.setZero(kNumFreqs, 1);
        chs_avg_speech_prob.setZero(kNumFreqs, 1);

        for (i = 0; i < kNumChannels; i++) {
            InputBufferRead(self->buf_rw_status[i].in_buffer.data(), &self->block_tf[kBlocks].in_td_block.col(i)[0],
                            kBlockLen, self->buf_rw_status[i].read_pos, self->buf_rw_status[i].write_pos,
                            self->buf_rw_status[i].rw_wrap, kBufferLen);
            MoveReadPositionBackward(-kBlockLapped, self->buf_rw_status[i].read_pos,
                                     self->buf_rw_status[i].write_pos, self->buf_rw_status[i].rw_wrap, kBufferLen);

            /* apply window to input block data */
            ApplyWindow1D(self->block_tf[kBlocks].in_td_block.data() + i * kBlockLen,
                          self->block_tf[kBlocks].in_td_block.data() + i * kBlockLen,
                          self->window.data(), kBlockLen);

            /* transform time domain data into frequency domain data ------ complex data */
            FftTransform(self->block_tf[kBlocks].in_td_block.data() + i * kBlockLen,
                         tmp_fd.data(), kBlockLen);

            /* save instant complex spectrum of each channel */
            tmp_fd_chs.col(i) = tmp_fd;
            tmp_fd_chs2.col(i) = tmp_fd;
            // compute power spectrum and angle spectrum
            ComputePowerDensitySp(self->block_tf[kBlocks].in_fd_block.col(i).data(), tmp_fd.data(), kNumFreqs);
            ComputeAngularSp(self->block_tf[kBlocks].ang_in_fd_block.col(i).data(), tmp_fd.data(), kNumFreqs);

            /* estimate snr, noise suppression gain, speech probability of each block for different channels */
            EstimateSnrAndGainEtc(&(self->block_tf[kBlocks].in_fd_block),
                                  self->snr_est, &self->noise_sps_est, &self->recur_avg_sps_est, &self->prior_snr,
                                  &self->post_snr, &self->gain, &self->speech_prob, self->frame_num, i);

            if (self->frame_num >= 5) {
                chs_avg_sf += self->recur_avg_sps_est.col(i);
                chs_avg_nf += self->noise_sps_est.col(i);
                chs_avg_speech_prob += self->speech_prob.col(i) / kNumChannels;
                tmp_fd_chs.col(i) = self->doa_core->sp_smoothing * tmp_fd_chs_prev.col(i) +
                                    (1 - self->doa_core->sp_smoothing) *
                                    tmp_fd_chs.col(i).cwiseProduct(self->gain.col(i));
                tmp_fd_chs_prev.col(i) = tmp_fd_chs.col(i);
            }

        } /* end of channels loop */

        /*
        voice detection:
        determine if current frame includes valid speech
        */
        VoiceDecisionWithinFrame(self->vadInst, chs_avg_sf.data(), chs_avg_nf.data(),
                                 chs_avg_speech_prob.data(), self->frame_num, self->vad_mode);

        self->is_speech = self->is_speech || self->vadInst->is_speech_frame;
        /*
        doa estimation and beamforming
        */
        if (self->is_speech) {
            /* doa estimation */
            DoaEstimationCore(&tmp_fd_chs, self->doa_core,
                              self->mic_array_info, self->doa_mode);
            self->doa = self->doa_core->doa;
            self->mic_array_info->azimuth = self->doa * deg_to_radian;
            printf("frame_num: %4d   doa:%.4f \n", self->frame_num, self->doa_core->doa);

            /* beamforming */
            // compute sound path for every scan_angle
            ComputeElemetsSoundPath(self->mic_array_info->location, self->mic_array_info->azimuth,
                                    self->mic_array_info->sound_path);
            BfCore(&tmp_fd_chs2, &tmp_fd_chs2, self->bf_core, self->mic_array_info, self->bf_mode);

        }


        /* post filter */
        tmp_fd_chs2 = tmp_fd_chs2.cwiseProduct(self->gain);
        PostFilterCore(&tmp_fd_chs2, &single_ch_out_sp, self->pf_core, self->cn_mode);


        /* transform frequency domain data into time domain data */
        InvFftTransform(single_ch_out_sp.data(), post_filter_output_block.data(), kBlockLen);

        /* apply window to output block data */
        ApplyWindow1D(post_filter_output_block.data(), post_filter_output_block.data(),
                      self->window.data(), kBlockLen);

        /* Add Frames --- Copy out_td_block into out_buffer */
        AddFrames(post_filter_output_buffer.data(), post_filter_output_block.data(),
                  first_frame_in_block, 0, kBlockLen);


        first_frame_in_block += kBlockShift;
        kBlocks++;
    } /* end of while loop */

    /* decide if current frame is within valid speech segment */
    ValidVoiceSegDecisionAlongFrames(self->vadInst, self->frame_num);
    self->is_speech_triggered = self->vadInst->is_speech_triggered;

    printf("frame_bumber:%4d   is_speech: %d  is_speech_triggered: %d \n",
           self->frame_num, self->is_speech, self->is_speech_triggered);


    // output foregoing kFrameLen samples of out_buffer as processed results
    // and update out_buffer for subsequent frames
    OutputSamplesAndUpdateBuffer(out_single_ch_data,
                                 post_filter_output_buffer.data(), kFrameLen, kBufferLen, kBlockInitialDelay);

    // Calculate new starting frames
    self->frame_offset = first_frame_in_block - kFrameLen;

    // increase sequence number of audio frame;
    self->frame_num++;
}

AMY_AUDIO_API void
ProcCoreMicArrayProcMultiOut(MicArrayProcInst *self, float **in_multi_chs_data, float **out_multi_chs_data) {
    int i = 0, j = 0, kBlocks = 0, first_frame_in_block = 0;

    /* read input audio data into audio ring buffer */
    for (i = 0; i < kNumChannels; i++) {
        InputBufferWrite(in_multi_chs_data[i], self->buf_rw_status[i].in_buffer.data(), kFrameLen,
                         self->buf_rw_status[i].read_pos, self->buf_rw_status[i].write_pos,
                         self->buf_rw_status[i].rw_wrap, kBufferLen);
    }

    /*
    based on input audio ring buffer, first_frame_in_block, etc,
    divide each frame data into blocks with kBlockLapped samples in adjacent blocks
    for each channel, process each frame block by block
    output data will be delayed by some samples
    */

    /* temporary variable to store complex spectrum */
    MatrixXcf tmp_fd, tmp_fd_chs;
    tmp_fd.setZero(kNumFreqs, 1);
    tmp_fd_chs.setZero(kNumFreqs, kNumChannels);


    first_frame_in_block = self->frame_offset;
    while (first_frame_in_block < kFrameLen) {
        for (i = 0; i < kNumChannels; i++) {
            InputBufferRead(self->buf_rw_status[i].in_buffer.data(), &self->block_tf[kBlocks].in_td_block.col(i)[0],
                            kBlockLen, self->buf_rw_status[i].read_pos, self->buf_rw_status[i].write_pos,
                            self->buf_rw_status[i].rw_wrap, kBufferLen);
            MoveReadPositionBackward(-kBlockLapped, self->buf_rw_status[i].read_pos,
                                     self->buf_rw_status[i].write_pos, self->buf_rw_status[i].rw_wrap, kBufferLen);

            /* apply window to input block data */
            ApplyWindow1D(self->block_tf[kBlocks].in_td_block.data() + i * kBlockLen,
                          self->block_tf[kBlocks].in_td_block.data() + i * kBlockLen,
                          self->window.data(), kBlockLen);

            /* transform time domain data into frequency domain data ------ complex data */
            FftTransform(self->block_tf[kBlocks].in_td_block.data() + i * kBlockLen,
                         tmp_fd.data(), kBlockLen);

            /* save instant complex spectrum of each channel */
            tmp_fd_chs.col(i) = tmp_fd;

            // compute power spectrum and angle spectrum
            ComputePowerDensitySp(self->block_tf[kBlocks].in_fd_block.col(i).data(), tmp_fd.data(), kNumFreqs);
            ComputeAngularSp(self->block_tf[kBlocks].ang_in_fd_block.col(i).data(), tmp_fd.data(), kNumFreqs);

            /* estimate snr, noise suppression gain, speech probability of each block for different channels */
            EstimateSnrAndGainEtc(&(self->block_tf[kBlocks].in_fd_block),
                                  self->snr_est, &self->noise_sps_est, &self->recur_avg_sps_est,
                                  &self->prior_snr, &self->post_snr, &self->gain, &self->speech_prob, self->frame_num,
                                  i);

            // apply Gainf to original audio signal spectrum
            ApplyGainToSp(self->gain.col(i).data(), self->block_tf[kBlocks].in_fd_block.col(i).data(),
                          self->block_tf[kBlocks].out_fd_block.col(i).data(), kNumFreqs);

            // construct complex spectrum
            ConstructComplexSpectrum(tmp_fd.data(), &self->block_tf[kBlocks].out_fd_block.col(i)[0],
                                     &self->block_tf[kBlocks].ang_in_fd_block.col(i)[0], kNumFreqs);

            // transform frequency domain data into time domain data
            InvFftTransform(tmp_fd.data(), &self->block_tf[kBlocks].out_td_block.col(i)[0], kBlockLen);

            // apply window to output block data
            ApplyWindow1D(self->block_tf[kBlocks].out_td_block.data() + i * kBlockLen,
                          self->block_tf[kBlocks].out_td_block.data() + i * kBlockLen,
                          self->window.data(), kBlockLen);


            // Add Frames --- Copy out_td_block into out_buffer
            AddFrames(self->buf_rw_status[i].out_buffer.data(),
                      &self->block_tf[kBlocks].out_td_block.col(i)[0],
                      first_frame_in_block, 0, kBlockLen);

        } /* end of channels loop */

        first_frame_in_block += kBlockShift;
        kBlocks++;
    } /* end of while loop */


    for (j = 0; j < kNumChannels; j++) {
        // output foregoing kFrameLen samples of out_buffer as processed results
        // and update out_buffer for subsequent frames
        OutputSamplesAndUpdateBuffer(out_multi_chs_data[j],
                                     self->buf_rw_status[j].out_buffer.data(), kFrameLen, kBufferLen,
                                     kBlockInitialDelay);
    }


    // Calculate new starting frames
    self->frame_offset = first_frame_in_block - kFrameLen;

    // increase sequence number of audio frame;
    self->frame_num++;

}


AMY_AUDIO_API void ComputePowerDensitySp(float *sp, const Complexf *xf, int len) {
    Complexf conj_temp(0.0f, 0.0f);
    for (int i = 0; i < len; i++) {
        conj_temp = Complexf(xf[i].real(), -xf[i].imag());
        sp[i] = (xf[i] * conj_temp).real();
    }
}

AMY_AUDIO_API void ComputeAngularSp(float *ang, const Complexf *xf, int len) {
    for (int i = 0; i < len; i++) {
        ang[i] = atan2f(xf[i].imag(), xf[i].real());
    }
}


AMY_AUDIO_API void
EstimateSnrAndGainEtc(MatrixXf *xf, VaNEstInst *psnr, MatrixXf *noise_est, MatrixXf *recur_avg_sp_est,
                      MatrixXf *prior_snr, MatrixXf *post_snr, MatrixXf *gain, MatrixXf *speech_prob, int frame_seq,
                      int ch_seq) {
    assert(ch_seq < kNumChannels && ch_seq >= 0);
    /* variable relative with snr estimation --- real-time update */
    /* local, will not be updated in frames, and will be reset every time */
    static MatrixXf InstantSmoothedPf = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf tmp_bkg_sp = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf tmp_Ef = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf min_Ef = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static ColVectorXf likelihood_ratio = ColVectorXf::Zero(kNumFreqs, 1);
    static ColVectorXf Ikl = ColVectorXf::Zero(kNumFreqs, 1);
    static MatrixXf SpeechPresenceProb = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf alpha_noise_sp_est = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->alpha_ns_prob_update);
    static MatrixXf RecurAvgPriorSNR = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static ColVectorXf P_Local = ColVectorXf::Zero(kNumFreqs, 1);
    static ColVectorXf P_Global = ColVectorXf::Zero(kNumFreqs, 1);
    float AvgRecurAvgPriorSNR = 0.0f;
    float P_Frame = 0.0f;
    float RecurAvgPriorSNR_Peak = 0.0f;
    static MatrixXf SpeechAbsenceProb = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static ColVectorXd vkPost_Prior = ColVectorXd::Zero(kNumFreqs, 1);
    static ColVectorXd expint_tmp = ColVectorXd::Zero(kNumFreqs, 1);
    static MatrixXf GainH1_Curr = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf GminMat = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->Gmin);


    /* global, will be updated in frames */
    static int blks = 0;
    static MatrixXf RecursiveEf_Prev = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf tmp_Ef_Prev = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf min_Ef_Prev = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf NoiseSpectrumEst_Prev = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf SpeechPresenceProb_Prev = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->prior_speech_presence);
    static MatrixXf PriorSNR_Prev = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->snr_prev_val);
    static MatrixXf RecurAvgPriorSNRPrev = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->snr_prev_val);
    static float RecurAvgPriorSNR_PeakPrev[kNumChannels] = {psnr->snr_prev_val, psnr->snr_prev_val};
    static float AvgRecurAvgPriorSNRPrev[kNumChannels] = {psnr->snr_prev_val, psnr->snr_prev_val};
    static MatrixXf GainH1_Prev = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->snr_prev_val);
    static MatrixXf PostSNR_Prev = MatrixXf::Constant(kNumFreqs, kNumChannels, psnr->snr_prev_val);




    /*
    leading NIS frames are used for estimating background noise
    accumulate power spectrum of leading NIS frames
    */
    if (frame_seq < NIS) {
        gain->setConstant(kNumFreqs, kNumChannels, psnr->Gmin);
        tmp_bkg_sp.col(ch_seq) += xf->col(ch_seq);
        blks++;
    }

    /* initial noise spectrum estimation */
    if (frame_seq == NIS) {
        MatrixXf smoothed_tmp_bkg_sp = MatrixXf::Zero(kNumFreqs, kNumChannels);
        /* average */
        tmp_bkg_sp.col(ch_seq) = tmp_bkg_sp.col(ch_seq).array() / blks * kNumChannels;
        /* smooth initial noise spectrum estimation */
        Smooth1Ddata(&tmp_bkg_sp.col(ch_seq)[0], &smoothed_tmp_bkg_sp.col(ch_seq)[0],
                     kNumFreqs, psnr->smoothed_win_len_of_sp,
                     (psnr->smooth_win_of_sp)->data());
        // initialize procedure variable for updating noise spectrum estimation
        RecursiveEf_Prev.col(ch_seq) = smoothed_tmp_bkg_sp.col(ch_seq);
        tmp_Ef_Prev.col(ch_seq) = smoothed_tmp_bkg_sp.col(ch_seq);
        min_Ef_Prev.col(ch_seq) = smoothed_tmp_bkg_sp.col(ch_seq);
        NoiseSpectrumEst_Prev.col(ch_seq) = tmp_bkg_sp.col(ch_seq);
    }

    /* real-time update process */
    if (frame_seq >= NIS) {
        // smooth original power spectrum for subsequent frames after foregoing NIS frames
        Smooth1Ddata(xf->col(ch_seq).data(), InstantSmoothedPf.col(ch_seq).data(),
                     kNumFreqs, psnr->smoothed_win_len_of_sp, psnr->smooth_win_of_sp->data());

        // update recursive value for RecursiveEf、min_Ef、tmp_Ef based on InstantSmoothedPf
        RecursiveAverage(RecursiveEf_Prev.col(ch_seq).data(), InstantSmoothedPf.col(ch_seq).data(),
                         recur_avg_sp_est->col(ch_seq).data(), kNumFreqs, psnr->alpha_sp);
        TranverseMinOf2Vecs(min_Ef_Prev.col(ch_seq).data(), recur_avg_sp_est->col(ch_seq).data(),
                            min_Ef.col(ch_seq).data(), kNumFreqs);
        TranverseMinOf2Vecs(tmp_Ef_Prev.col(ch_seq).data(), recur_avg_sp_est->col(ch_seq).data(),
                            tmp_Ef.col(ch_seq).data(), kNumFreqs);

        if (!((frame_seq + 1) % psnr->noise_minsp_est_win_len)) {
            TranverseMinOf2Vecs(tmp_Ef_Prev.col(ch_seq).data(), recur_avg_sp_est->col(ch_seq).data(),
                                min_Ef.col(ch_seq).data(), kNumFreqs);
            tmp_Ef.col(ch_seq) = recur_avg_sp_est->col(ch_seq);
        }

        // compute likelihood ratio for frequency-wise power spectrum
        ComputeDotQuotientOf2Vecs(recur_avg_sp_est->col(ch_seq).data(), min_Ef.col(ch_seq).data(),
                                  likelihood_ratio.data(), kNumFreqs);

        // compute frequency-wise probability of existence of speech for current frame
        CompareVecWithScalar(likelihood_ratio.data(), psnr->likelihood_ratio_threshold, Ikl.data(), kNumFreqs);
        RecursiveAverage(SpeechPresenceProb_Prev.col(ch_seq).data(), Ikl.data(),
                         SpeechPresenceProb.col(ch_seq).data(), kNumFreqs, psnr->alpha_speech_presence_prob);

        // update smoothing factor for updating noise spectrum estimation
        // to be considered deeply
        RecurAveVecWithScalar(SpeechPresenceProb.col(ch_seq).data(), psnr->alpha_ns_prob_update,
                              alpha_noise_sp_est.col(ch_seq).data(), kNumFreqs, psnr->alpha_ns_prob_update);
        RecursiveAverageWithSmoothVec(NoiseSpectrumEst_Prev.col(ch_seq).data(), InstantSmoothedPf.col(ch_seq).data(),
                                      noise_est->col(ch_seq).data(),
                                      kNumFreqs, alpha_noise_sp_est.col(ch_seq).data());

        // update relative variable for noise spectrum estimation
        CopyValsOfVecsToOther(recur_avg_sp_est->col(ch_seq).data(), RecursiveEf_Prev.col(ch_seq).data(), kNumFreqs);
        CopyValsOfVecsToOther(min_Ef.col(ch_seq).data(), min_Ef_Prev.col(ch_seq).data(), kNumFreqs);
        CopyValsOfVecsToOther(tmp_Ef.col(ch_seq).data(), tmp_Ef_Prev.col(ch_seq).data(), kNumFreqs);
        CopyValsOfVecsToOther(SpeechPresenceProb.col(ch_seq).data(), SpeechPresenceProb_Prev.col(ch_seq).data(),
                              kNumFreqs);
        CopyValsOfVecsToOther(noise_est->col(ch_seq).data(), NoiseSpectrumEst_Prev.col(ch_seq).data(), kNumFreqs);

        /* speech estimation */
        // compute recursive average of prior snr
        RecursiveAverage(RecurAvgPriorSNRPrev.col(ch_seq).data(), PriorSNR_Prev.col(ch_seq).data(),
                         RecurAvgPriorSNR.col(ch_seq).data(), kNumFreqs, psnr->beta_prior_snr);

        // smooth recursive average of prior snr locally
        Smooth1Ddata(RecurAvgPriorSNR.col(ch_seq).data(), P_Local.data(),
                     kNumFreqs, psnr->smoothed_local_win_len, psnr->smoothed_local_win->data());

        // smooth recursive average of prior snr globally
        Smooth1Ddata(RecurAvgPriorSNR.data() + ch_seq * kNumFreqs, P_Global.data(),
                     kNumFreqs, psnr->smoothed_global_win_len, psnr->smoothed_global_win->data());

        // constrain local and global snr
        ConstrainLocalGlobalSnr(P_Local.data(), P_Global.data(), psnr->sigma_min, psnr->sigma_max, kNumFreqs);

        // compute average snr of leading frequency bins of certain number
        AvgRecurAvgPriorSNR = ComputeVecHeadAvg(RecurAvgPriorSNR.col(ch_seq).data(), 0, (kNumFreqs - 1) / 2 + 1);

        // Compute snr feature for whole frame and update snr peak value
        ComputePFrameAndSnrPeak(&P_Frame, &RecurAvgPriorSNR_Peak, RecurAvgPriorSNR_PeakPrev, AvgRecurAvgPriorSNR,
                                AvgRecurAvgPriorSNRPrev, psnr->sigma_min, psnr->sigma_max, psnr->sigma_pmin,
                                psnr->sigma_pmax, ch_seq);

        // compute absence probability of speech
        ComputeSpeechAbsenceProb(SpeechAbsenceProb.col(ch_seq).data(), P_Local.data(),
                                 P_Global.data(), P_Frame, kNumFreqs);

        // Limit maximum absence probability of speech to 0.95
        GetMinOfVecAndScalar(SpeechAbsenceProb.col(ch_seq).data(), SpeechAbsenceProb.col(ch_seq).data(),
                             psnr->max_speech_absence_prob, kNumFreqs);

        // update relative variable for speech probability
        AvgRecurAvgPriorSNRPrev[ch_seq] = AvgRecurAvgPriorSNR;
        CopyValsOfVecsToOther(RecurAvgPriorSNR.col(ch_seq).data(), RecurAvgPriorSNRPrev.col(ch_seq).data(), kNumFreqs);

        // compute posterior snr
        ComputeDotQuotientOf2Vecs(xf->col(ch_seq).data(), noise_est->col(ch_seq).data(),
                                  post_snr->col(ch_seq).data(), kNumFreqs);

        // estimate and update gain, prior and posterior snr
        ComputePriorSNRCurr(prior_snr->col(ch_seq).data(), GainH1_Prev.col(ch_seq).data(),
                            PostSNR_Prev.col(ch_seq).data(), post_snr->col(ch_seq).data(), psnr->compromising,
                            kNumFreqs);
        ComputeVkPost_Prior(vkPost_Prior.data(), prior_snr->col(ch_seq).data(),
                            post_snr->col(ch_seq).data(), kNumFreqs);

        Exponential_Integral_Ei_vector(vkPost_Prior.data(), expint_tmp.data(), kNumFreqs);

        ComputeGain_Curr(GainH1_Curr.col(ch_seq).data(), prior_snr->col(ch_seq).data(),
                         expint_tmp.data(), kNumFreqs);

        // update Gain and SNR estiamtions
        CopyValsOfVecsToOther(prior_snr->col(ch_seq).data(), PriorSNR_Prev.col(ch_seq).data(), kNumFreqs);
        CopyValsOfVecsToOther(post_snr->col(ch_seq).data(), PostSNR_Prev.col(ch_seq).data(), kNumFreqs);
        GetMaxOfVecAndScalar(GainH1_Curr.col(ch_seq).data(), GainH1_Prev.col(ch_seq).data(), psnr->Gmin, kNumFreqs);

        // compute speech presence probability of current frame
        ComputeSpeechPresenceProb(SpeechAbsenceProb.col(ch_seq).data(), prior_snr->col(ch_seq).data(),
                                  vkPost_Prior.data(), speech_prob->col(ch_seq).data(), kNumFreqs);

        // compute final frequency-wise gain
        ComputeFinalGain(gain->col(ch_seq).data(), GainH1_Curr.col(ch_seq).data(),
                         speech_prob->col(ch_seq).data(), psnr->Gmin, kNumFreqs);

    } /* if (frame_seq >= NIS) */



}

AMY_AUDIO_API void RecursiveAverage(float *prev, float *curr, float *dst, int len, float factor) {
    for (int i = 0; i < len; i++) {
        dst[i] = factor * prev[i] + (1 - factor) * curr[i];
    }
}

AMY_AUDIO_API void TranverseMinOf2Vecs(float *vec1, float *vec2, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] < vec2[i] ? vec1[i] : vec2[i];
    }
}

AMY_AUDIO_API void TranverseMaxOf2Vecs(float *vec1, float *vec2, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] > vec2[i] ? vec1[i] : vec2[i];
    }
}

AMY_AUDIO_API void ComputeDotQuotientOf2Vecs(float *num_vec, float *den_vec, float *dst, int len) {
    float denominator = 0.0f;
    for (int i = 0; i < len; i++) {
        denominator = den_vec[i] > 1e-6f ? den_vec[i] : 1e-6f;
        dst[i] = num_vec[i] / denominator;
    }
}

AMY_AUDIO_API void ComputeDotProductOf2Vecs(float *vec1, float *vec2, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] * vec2[i];
    }
}

AMY_AUDIO_API void CompareVecWithScalar(float *vec1, float scalar, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] > scalar ? 1.0f : 0.0f;
    }
}

AMY_AUDIO_API void RecurAveVecWithScalar(float *vec, float scalar, float *dst, int len, float factor) {
    for (int i = 0; i < len; i++) {
        dst[i] = scalar + (1 - factor) * vec[i];
    }
}

AMY_AUDIO_API void RecursiveAverageWithSmoothVec(float *prev, float *curr, float *dst, int len, float *factor) {
    for (int i = 0; i < len; i++) {
        dst[i] = factor[i] * prev[i] + (1 - factor[i]) * curr[i];
    }
}

AMY_AUDIO_API void CopyValsOfVecsToOther(float *src, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = src[i];
    }
}

AMY_AUDIO_API void ConstrainLocalGlobalSnr(float *local_snr, float *global_snr, float min_th, float max_th, int len) {
    for (int i = 0; i < kNumFreqs; i++) {
        if (local_snr[i] <= min_th)
            local_snr[i] = 0.0f;
        else if (local_snr[i] >= max_th)
            local_snr[i] = 1.0f;
        else
            local_snr[i] = (log10f(local_snr[i] / min_th)) /
                           (log10f(max_th / min_th));

        if (global_snr[i] <= min_th)
            global_snr[i] = 0;
        else if (global_snr[i] >= max_th)
            global_snr[i] = 1;
        else
            global_snr[i] = (log10f(global_snr[i] / min_th)) /
                            (log10f(max_th / min_th));
    }
}

AMY_AUDIO_API float ComputeVecHeadAvg(float *vec, int start_ind, int end_ind) {
    float tmp = 0.0f;
    for (int i = start_ind; i < end_ind; i++) {
        tmp += vec[i];
    }
    return tmp / (end_ind - start_ind);
}

AMY_AUDIO_API void ComputePFrameAndSnrPeak(float *p_frame, float *snr_peak, float *snr_peak_prev, float seg_avg_snr,
                                           float *seg_avg_snr_prev, float snr_min, float snr_max, float snr_pmin,
                                           float snr_pmax, int ch_seq) {
    float th1 = snr_min * powf(10, snr_peak_prev[ch_seq] / 20);
    float th2 = snr_max * powf(10, snr_peak_prev[ch_seq] / 20);
    float mu_l = 0.0f;
    if (seg_avg_snr <= th1)
        mu_l = 0;
    else if (seg_avg_snr >= th2)
        mu_l = 1;
    else
        mu_l = log10f(seg_avg_snr / th1) / log10f(snr_max / snr_min);

    // compute P_Frame
    if (seg_avg_snr <= snr_min)
        *p_frame = 0.0f;
    else if (seg_avg_snr <= seg_avg_snr_prev[ch_seq])
        *p_frame = mu_l;
    else {
        *snr_peak = seg_avg_snr >= snr_pmin ? seg_avg_snr : snr_pmin;
        *snr_peak = *snr_peak > snr_pmax ? snr_pmax : *snr_peak;
        snr_peak_prev[ch_seq] = *snr_peak;
        *p_frame = 1;
    }
}


AMY_AUDIO_API void ComputeSpeechAbsenceProb(float *absence_prob, float *local_snr,
                                            float *global_snr, float frame_snr, int len) {
    for (int i = 0; i < len; i++) {
        absence_prob[i] = 1 - local_snr[i] * global_snr[i] * frame_snr;
    }
}


AMY_AUDIO_API void GetMinOfVecAndScalar(float *vec, float *min_vec, float scalar, int len) {
    for (int i = 0; i < len; i++) {
        min_vec[i] = vec[i] < scalar ? vec[i] : scalar;
    }
}

AMY_AUDIO_API void GetMaxOfVecAndScalar(float *vec, float *max_vec, float scalar, int len) {
    for (int i = 0; i < len; i++) {
        max_vec[i] = vec[i] > scalar ? vec[i] : scalar;
    }
}


AMY_AUDIO_API void ComputePriorSNRCurr(float *prior_snr, float *gain_prev,
                                       float *post_snr_prev, float *post_snr_curr, float compromising, int len) {
    for (int i = 0; i < len; i++) {
        prior_snr[i] = compromising * powf(gain_prev[i], 2) * post_snr_prev[i] +
                       (1 - compromising) * std::max(post_snr_curr[i] - 1, 0.0f);
    }
}

AMY_AUDIO_API void ComputeVkPost_Prior(double *vk, float *prior_snr_curr, float *post_snr_curr, int len) {
    for (int i = 0; i < len; i++) {
        vk[i] = -1.0 * prior_snr_curr[i] * post_snr_curr[i] / (1 + prior_snr_curr[i]);
    }
}

AMY_AUDIO_API void ComputeGain_Curr(float *gain_curr, float *prior_snr_curr, double *exp_integeral, int len) {
    float exp_tmp = 0.0f;
    for (int i = 0; i < len; i++) {
        exp_tmp = float(std::exp(-1.0f * exp_integeral[i] / 2));
        gain_curr[i] = prior_snr_curr[i] / (1 + prior_snr_curr[i]) * exp_tmp;
    }
}

AMY_AUDIO_API void ComputeSpeechPresenceProb(float *speech_absence_prob, float *prior_snr_curr,
                                             double *vk_post_prior, float *speech_presence_prob, int len) {
    for (int i = 0; i < len; i++) {
        speech_presence_prob[i] = 1 + speech_absence_prob[i] / (1 - speech_absence_prob[i]) *
                                      (1 + prior_snr_curr[i]) * float(std::exp(vk_post_prior[i]));
        speech_presence_prob[i] = 1 / speech_presence_prob[i];
    }

}

AMY_AUDIO_API void ComputeFinalGain(float *final_gain, float *gain_curr,
                                    float *speech_presence_prob, float gain_min, int len) {
    for (int i = 0; i < len; i++) {
        final_gain[i] = powf(gain_curr[i], speech_presence_prob[i]) *
                        powf(gain_min, 1 - speech_presence_prob[i]);
    }
}

AMY_AUDIO_API void ApplyGainToSp(float *gain, float *original_sp, float *processed_sp, int len) {
    for (int i = 0; i < len; i++) {
        processed_sp[i] = original_sp[i] * powf(gain[i], 2);
    }
}


AMY_AUDIO_API void VoiceDecisionWithinFrame(VoiceDetInst *vad, float *avg_sf, float *avg_nf,
                                            float *voice_prob, int frame_seq, int vad_mode) {
    /* fundamental variables */
    static float bkg_noise_energy = 0.0f;
    static float bkg_noise_energy_entropy_ratio = 0.0f;
    static float instant_energy = 0.0f;
    static float instant_energy_entropy_ratio = 0.0f;
    ColVectorXf smoothed_voice_prob = ColVectorXf::Zero(kNumFreqs, 1);
    static ColVectorXf smoothed_voice_prob_prev = ColVectorXf::Constant(kNumFreqs, 1, vad->initial_voice_prob);
    static float avg_voice_prob_largers = 0.0f;
    static float avg_voice_prob_largers_prev = 0.0f;

    if (frame_seq < NIS) {
        vad->is_speech_frame = false;
        vad->is_speech_triggered = false;
    } else {
        switch (vad_mode) {
            case 0:
                /* process speech probability */
                Smooth1Ddata(voice_prob, smoothed_voice_prob.data(),
                             kNumFreqs, vad->smoothed_win_len_of_sp, vad->smooth_win_of_sp->data());
                RecursiveAverage(smoothed_voice_prob_prev.data(), smoothed_voice_prob.data(),
                                 smoothed_voice_prob.data(), kNumFreqs, vad->smoothing);
                CopyValsOfVecsToOther(smoothed_voice_prob.data(), smoothed_voice_prob_prev.data(), kNumFreqs);
                /* sort smoothed_voice_prob with descending order */
                qsort(smoothed_voice_prob.data(), kNumFreqs, sizeof(float), cmp);
                /* compute mean value of leading LargerNum larger values*/
                avg_voice_prob_largers = ComputeVecHeadAvg(smoothed_voice_prob.data(), 0, vad->LargerNum);
                avg_voice_prob_largers = vad->smoothing * avg_voice_prob_largers_prev +
                                         (1 - vad->smoothing) * avg_voice_prob_largers;
                avg_voice_prob_largers_prev = avg_voice_prob_largers;

                /* vad decision */
                vad->is_speech_frame = avg_voice_prob_largers > vad->voice_prob_threshold ?
                                       true : false;

                break;
            case 1:
                /* compute average noise energy within certain frequency range */
                bkg_noise_energy = ComputeIntervalAverage(avg_nf, vad->frequency_start_ind, vad->frequency_end_ind);
                /* compute average constant energy with certain frequency range for current frame */
                instant_energy = ComputeIntervalAverage(avg_sf, vad->frequency_start_ind, vad->frequency_end_ind);
                vad->is_speech_frame = instant_energy > bkg_noise_energy * vad->beyond_ratio_energy ?
                                       true : false;
                break;
            case 2:
                /* compute noise energy entropy ratio with certain frequency range */
                bkg_noise_energy_entropy_ratio = ComputeLogEnergyEntropyRatio(avg_nf, vad->frequency_start_ind,
                                                                              vad->frequency_end_ind,
                                                                              vad->alpha_log_sp);
                /* compute energy entropy ratio with certain frequency range for current frame */
                instant_energy_entropy_ratio = ComputeLogEnergyEntropyRatio(avg_sf, vad->frequency_start_ind,
                                                                            vad->frequency_end_ind, vad->alpha_log_sp);
                vad->is_speech_frame =
                        instant_energy_entropy_ratio > bkg_noise_energy_entropy_ratio * vad->beyond_ratio_ee ?
                        true : false;
                break;
            default:
                printf("vad_mode value must be 0, 1, or 2!!!\n");
        } /* end of switch */
    } /* end of if condition */
}


AMY_AUDIO_API void ValidVoiceSegDecisionAlongFrames(VoiceDetInst *vad, int frame_seq) {
    int non_speech_to_speech_ret = 0;
    int speech_to_non_speech_ret = 0;
    static MatrixXi non_speech_to_speech_list = MatrixXi::Zero(vad->non_speech_to_speech_frame_len, 1);
    static MatrixXi speech_to_non_speech_list = MatrixXi::Ones(vad->speech_to_non_speech_frame_len, 1);
    // decide whether current frame is in valid speech segments or not, namely is_speech_triggered is true or not
    int valid_nums = 0, sum_true = 0;
    if (!vad->is_speech_triggered) {
        /* judge length valid decision window */
        valid_nums =
                frame_seq < vad->non_speech_to_speech_frame_len ? frame_seq + 1 : vad->non_speech_to_speech_frame_len;
        /* judge element to be updated in valid decision window */
        non_speech_to_speech_ret = frame_seq % vad->non_speech_to_speech_frame_len;
        /* update value of certain element in valid decision window */
        non_speech_to_speech_list(non_speech_to_speech_ret) = vad->is_speech_frame;
        /* compute number of speech frame in valid decision window */
        sum_true = ComputeSumOfBoolVec(non_speech_to_speech_list.data(), valid_nums);

        /* start of valid speech segment decision */
        if (1.0f * sum_true / vad->non_speech_to_speech_frame_len > vad->speech_start_threshold) {
            vad->is_speech_triggered = true;
            non_speech_to_speech_list.setZero();
        }
    } else {
        valid_nums =
                frame_seq < vad->speech_to_non_speech_frame_len ? frame_seq + 1 : vad->speech_to_non_speech_frame_len;
        speech_to_non_speech_ret = frame_seq % vad->speech_to_non_speech_frame_len;
        speech_to_non_speech_list(speech_to_non_speech_ret) = vad->is_speech_frame;
        sum_true = ComputeSumOfBoolVec(speech_to_non_speech_list.data(), valid_nums);

        /* end of valid speech segment decision */
        if (1.0f * sum_true / vad->speech_to_non_speech_frame_len < vad->speech_finished_threshold) {
            vad->is_speech_triggered = false;
            speech_to_non_speech_list.setOnes();
        }
    }
}


AMY_AUDIO_API int ComputeSumOfBoolVec(int *bvec, int len) {
    int tmp = 0;
    for (int i = 0; i < len; i++) {
        tmp += bvec[i];
    }
    return tmp;
}

AMY_AUDIO_API void DoaEstimationCore(MatrixXcf *chs_com_sp, DoaInst *pdoa, ArrayProps *mic_info, int doa_mode) {
    /*
    use complex spectrum to perform beam scan
    estimate direction of arrival of sound source of speaker
    find maximum value corresponding to spatial spectrum as doa estimation
    */

    float frequency_scan = 0.0f;   /* instant sweep frequency */

    static int scan_angle_num = mic_info->is_linear ? (int((pdoa->coarse_ang_max1 - pdoa->coarse_ang_start) /
                                                           pdoa->coarse_ang_span) + 1) : (
                                        int((pdoa->coarse_ang_max2 - pdoa->coarse_ang_start) /
                                            pdoa->coarse_ang_span) + 1);


    MatrixXcf array_steer_vecter, freq_wise_complex_spectrum;
    array_steer_vecter.setZero(kNumChannels, 1);
    freq_wise_complex_spectrum.setZero(kNumChannels, 1);

    MatrixXf Pf_out = MatrixXf::Zero(scan_angle_num, 1);
    Complexf RSTCM[kNumChannels][kNumChannels];
    pdoa->doa = 0.0f;
    int k1 = 0, k2 = 0, max_index = 0;

    /* sweep angle */
    for (k1 = 0; k1 < scan_angle_num; k1++) {
        SetComplexfVec2Zeros(RSTCM, kNumChannels);
        mic_info->azimuth = (pdoa->coarse_ang_start + k1 * pdoa->coarse_ang_span)
                            * deg_to_radian;
        // compute sound path for every scan_angle
        ComputeElemetsSoundPath(mic_info->location, mic_info->azimuth,
                                mic_info->sound_path);

        /* sweep frequency */
        for (k2 = pdoa->frequency_start_ind; k2 <= pdoa->frequency_end_ind; k2++) {
            frequency_scan = pdoa->freq_bin[k2];
            // compute array steering vector based on sound_path and frequency
            DoaEstComputeSteerVec(mic_info->sound_path, array_steer_vecter.data(),
                                  frequency_scan);
            // extract complex spectrum of certain frequency
            freq_wise_complex_spectrum.col(0) = chs_com_sp->row(k2);

            switch (doa_mode) {
                /* conventional beamforming */
                case 0:
                    Pf_out(k1) += ComputeOutSpatialSpectrumCbf(array_steer_vecter.data(),
                                                               freq_wise_complex_spectrum.data(), kNumChannels);
                    break;
                    /* minimum variance distortionless response beamforming */
                case 1:
                    Pf_out(k1) += ComputeOutSpatialSpectrumMvdr(array_steer_vecter.data(),
                                                                freq_wise_complex_spectrum.data(), kNumChannels,
                                                                pdoa->loading_factor);
                    break;
                    /* steered minimum variance beamforming */
                case 2:
                    ComputeStmvNormCov(array_steer_vecter.data(), freq_wise_complex_spectrum.data(), RSTCM);
                    break;
                default:
                    printf("doa_mode value must be 0, 1, or 2!!!\n");
            } /* end of switch */
            if (doa_mode == 2)
                Pf_out(k1) = ComputeOutSpatialSpectrumStmv(RSTCM, kNumChannels, pdoa->loading_factor);
        }
    }

    /* find doa estimate corresponding to maximum spatial spectrum */
    Pf_out.col(0).maxCoeff(&max_index);
    pdoa->doa = pdoa->coarse_ang_start + max_index * pdoa->coarse_ang_span;

}

AMY_AUDIO_API void SetComplexfVec2Zeros(Complexf(*vec2D)[kNumChannels], int rows) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < kNumChannels; j++)
            vec2D[i][j] = Complexf(0.0f, 0.0f);
}


AMY_AUDIO_API void DoaEstComputeSteerVec(float *sound_path, Complexf *steer_vcetor, float frequency) {
    float tao = 0.0f, real = 0.0f, imag = 0.0f;
    for (int i = 0; i < ArrayElements; i++) {
        tao = sound_path[i] / sound_speed;
        real = cosf(-2 * M_PI_F * frequency * tao);
        imag = sinf(-2 * M_PI_F * frequency * tao);
        steer_vcetor[i] = Complexf(real, imag);
    }
}

AMY_AUDIO_API Complexf ComputeCompVecConjDotProduct(const Complexf *vec1, const Complexf *vec2, int len) {
    Complexf ctmp(0.0f, 0.0f), conj_steer_vector(0.0f, 0.0f);
    for (int i = 0; i < len; i++) {
        conj_steer_vector = Complexf(vec1[i].real(), -vec1[i].imag());
        ctmp += conj_steer_vector * vec2[i];
    }
    return ctmp;
}

AMY_AUDIO_API float
ComputeOutSpatialSpectrumCbf(const Complexf *steer_vcetor, const Complexf *In_FreqSpectrum, int len) {
    float tmp = 0.0f;
    Complexf ctmp(0.0f, 0.0f), conj_steer_vector(0.0f, 0.0f);
    ctmp = ComputeCompVecConjDotProduct(steer_vcetor, In_FreqSpectrum, len);
    tmp = ctmp.real() * ctmp.real() + ctmp.imag() * ctmp.imag();
    return tmp;
}


AMY_AUDIO_API float NormComplex(Complexf &a) {
    return sqrtf(a.real() * a.real() + a.imag() * a.imag());
}


AMY_AUDIO_API Complexf ComputeCompVecDotProduct(const Complexf *vec1, const Complexf *vec2, int len) {
    Complexf ctmp(0.0f, 0.0f);
    for (int i = 0; i < len; i++) {
        ctmp += vec1[i] * vec2[i];
    }
    return ctmp;
}


AMY_AUDIO_API float
ComputeOutSpatialSpectrumMvdr(const Complexf *steer_vcetor, const Complexf *In_FreqSpectrum, int len,
                              float load_factor) {
    float tmp = 0.0f;
    Complexf ctmp(0.0f, 0.0f), conj_steer_vector(0.0f, 0.0f);
    MatrixXcf Rxf = MatrixXcf::Zero(kNumChannels, kNumChannels);
    MatrixXcf invRxf = MatrixXcf::Zero(kNumChannels, kNumChannels);
    int i = 0, j = 0;
    /* compute spatial spectrum covariance matrix with loading factor */
    for (i = 0; i < kNumChannels; i++)
        for (j = 0; j < kNumChannels; j++) {
            conj_steer_vector = Complexf(In_FreqSpectrum[j].real(), -In_FreqSpectrum[j].imag());
            if (i == j)
                Rxf(i, j) = In_FreqSpectrum[i] * conj_steer_vector + Complexf(load_factor, 0.0f);
            else
                Rxf(i, j) = In_FreqSpectrum[i] * conj_steer_vector;
        }
    invRxf = Rxf.inverse();
    /* compute output spatial spectrum */
    Complexf ctmp_mid[kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
        ctmp_mid[i] = ComputeCompVecConjDotProduct(steer_vcetor, &invRxf(0, i), kNumChannels);

    ctmp = ComputeCompVecDotProduct(ctmp_mid, steer_vcetor, kNumChannels);
    tmp = 1.0f / NormComplex(ctmp);

    return tmp;
}


AMY_AUDIO_API void ComputeStmvNormCov(const Complexf *vec1, const Complexf *vec2, Complexf(*norm_mat)[kNumChannels]) {
    Complexf ctmp_mid[kNumChannels], conj_tmp(0.0f, 0.0f);
    for (int i = 0; i < kNumChannels; i++) {
        conj_tmp = Complexf(vec1[i].real(), -vec1[i].imag());
        ctmp_mid[i] = conj_tmp * vec2[i];
    }

    Complexf tmp_mat[kNumChannels][kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
        for (int j = 0; j < kNumChannels; j++) {
            conj_tmp = Complexf(vec2[j].real(), -vec2[j].imag());
            tmp_mat[i][j] = ctmp_mid[i] * conj_tmp;
        }

    Complexf tmp_mat2[kNumChannels][kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
        for (int j = 0; j < kNumChannels; j++) {
            tmp_mat2[i][j] = tmp_mat[i][j] * vec1[j];
        }

    for (int i = 0; i < kNumChannels; i++)
        for (int j = 0; j < kNumChannels; j++)
            norm_mat[i][j] += tmp_mat2[i][j];


}


AMY_AUDIO_API float ComputeOutSpatialSpectrumStmv(Complexf(*norm_mat)[kNumChannels], int len, float load_factor) {
    float tmp = 0.0f;
    Complexf ctmp(0.0f, 0.0f);
    ColVectorXcf OnesVec = ColVectorXcf::Ones(len);
    assert(kNumChannels == len);
    MatrixXcf cov_mat = MatrixXcf::Zero(kNumChannels, kNumChannels);
    MatrixXcf inv_cov_mat = MatrixXcf::Zero(kNumChannels, kNumChannels);
    for (int i = 0; i < kNumChannels; i++) {
        for (int j = 0; j < kNumChannels; j++) {
            norm_mat[i][j] = i == j ? norm_mat[i][j] + load_factor : norm_mat[i][j];
            cov_mat(i, j) = norm_mat[i][j];
        }
    }

    inv_cov_mat = cov_mat.inverse();

    /* compute output spatial spectrum */
    Complexf ctmp_mid[kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
        ctmp_mid[i] = ComputeCompVecConjDotProduct(OnesVec.data(), &inv_cov_mat(0, i), kNumChannels);

    ctmp = ComputeCompVecDotProduct(ctmp_mid, OnesVec.data(), kNumChannels);

    tmp = 1.0f / NormComplex(ctmp);

    return tmp;
}

AMY_AUDIO_API void BfComputeSteerVec(float *sound_path, Complexf *steer_vcetor, float frequency) {
    float tao = 0.0f, real = 0.0f, imag = 0.0f;
    for (int i = 0; i < ArrayElements; i++) {
        tao = (sound_path[i] - sound_path[0]) / sound_speed;
        real = cosf(-2 * M_PI_F * frequency * tao);
        imag = sinf(-2 * M_PI_F * frequency * tao);
        steer_vcetor[i] = Complexf(real, imag);
    }
}


AMY_AUDIO_API void
BfCore(MatrixXcf *in_chs_com_sp, MatrixXcf *out_chs_com_sp, bfInst *pbf, ArrayProps *mic_info, int bf_mode) {
    MatrixXcf array_steer_vecter, freq_wise_complex_spectrum; /* to store complex spectrum */
    array_steer_vecter.setZero(kNumChannels, 1);
    freq_wise_complex_spectrum.setZero(kNumChannels, 1);
    float frequency_scan = 0.0f;
    int k = 0;
    for (k = 0; k < kNumFreqs; k++) {
        frequency_scan = pbf->freq_bin[k];
        // compute array steering vector
        BfComputeSteerVec(mic_info->sound_path, array_steer_vecter.data(),
                          frequency_scan);
        // extract complex spectrum of certain frequency
        freq_wise_complex_spectrum = in_chs_com_sp->row(k).transpose();
        // conventional beamforming
        BfPhaseCompensation(array_steer_vecter.data(), freq_wise_complex_spectrum.data(),
                            freq_wise_complex_spectrum.data(), kNumChannels);
        out_chs_com_sp->row(k) = freq_wise_complex_spectrum.transpose();
    }
}


AMY_AUDIO_API void BfPhaseCompensation(const Complexf *steer_vec, Complexf *src_sp, Complexf *dst_sp, int len) {
    Complexf conj_tmp(0.0f, 0.0f);
    for (int i = 0; i < len; i++) {
        conj_tmp = Complexf(steer_vec[i].real(), -steer_vec[i].imag());
        dst_sp[i] = conj_tmp * src_sp[i];
    }
}


AMY_AUDIO_API void PostFilterCore(MatrixXcf *in_chs_com_sp, MatrixXcf *single_out_ch_com_sp, PfInst *ppf, int cn_mode) {
    static MatrixXf auto_sp_prev = MatrixXf::Zero(kNumFreqs, kNumChannels);
    MatrixXf auto_sp = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf inter_sp_prev = MatrixXf::Zero(kNumFreqs, ppf->Num_Inters);
    MatrixXf inter_sp = MatrixXf::Zero(kNumFreqs, ppf->Num_Inters);
    MatrixXf speech_sp = MatrixXf::Zero(kNumFreqs, 1);
    MatrixXf speech_plus_noise_sp = MatrixXf::Zero(kNumFreqs, 1);
    MatrixXf post_filter_coeff = MatrixXf::Zero(kNumFreqs, 1);


    /* compute self power spectrum of each channel */
    for (int kk = 0; kk < kNumChannels; kk++) {
        ComputePowerDensity(auto_sp_prev.col(kk).data(), auto_sp.col(kk).data(),
                            in_chs_com_sp->col(kk).data(), in_chs_com_sp->col(kk).data(), ppf->smooth_alpha, kNumFreqs);
    }

    /* compute cross power spectrum of inter channels */
    int k = 0, k1 = 0, k2 = 0;
    for (k1 = 0; k1 < kNumChannels - 1; k1++)
        for (k2 = k1 + 1; k2 < kNumChannels; k2++) {
            ComputePowerDensity(inter_sp_prev.col(k).data(), inter_sp.col(k).data(),
                                in_chs_com_sp->col(k1).data(), in_chs_com_sp->col(k2).data(), ppf->smooth_alpha,
                                kNumFreqs);

            k++;
        }

    /* estimate speech power spectrum and speech plus noise power spectrum */
    ComputeAveragePowerDensity(auto_sp.data(), speech_plus_noise_sp.data(), kNumChannels, kNumFreqs);
    ComputeAveragePowerDensity(inter_sp.data(), speech_sp.data(), ppf->Num_Inters, kNumFreqs);

    /* compute post filter coeffcient */
    ComputePostFilterCoeff(speech_sp.data(), speech_plus_noise_sp.data(), post_filter_coeff.data(), kNumFreqs);
    if (cn_mode == 0) {
        GetCwiseMaxValOfVecAndScalar(post_filter_coeff.data(), ppf->min_postfilter_coeff, post_filter_coeff.data(),
                                     kNumFreqs);
    } else if (cn_mode == 1) {
        GetCwiseMaxValOf2Vecs(post_filter_coeff.data(), ppf->min_postfilter_coeff_cn->data(), post_filter_coeff.data(),
                              kNumFreqs);
    }

    GetCwiseMinValOfVecAndScalar(post_filter_coeff.data(), 1.0f, post_filter_coeff.data(), kNumFreqs);

    /* apply postfilter to get output complex spectrum: first way */
//	ComputeAverageComSpectrum(in_chs_com_sp->data(), single_out_ch_com_sp->data(), kNumChannels, kNumFreqs);
//	ApplyPostFilterCoeff(single_out_ch_com_sp->data(), single_out_ch_com_sp->data(), post_filter_coeff.data(), kNumFreqs);

    /* apply postfilter to get output complex spectrum: second way */
    static MatrixXf square_sp = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf angle_sp = MatrixXf::Zero(kNumFreqs, kNumChannels);
    static MatrixXf avg_square_sp = MatrixXf::Zero(kNumFreqs, 1);
    square_sp = in_chs_com_sp->array().abs().pow(2);
    angle_sp = in_chs_com_sp->array().arg();
    avg_square_sp = square_sp.rowwise().mean();

    ConstructComSpFromEsp_Psp(single_out_ch_com_sp->data(), avg_square_sp.col(0).data(), angle_sp.data(), kNumFreqs);
    ApplyPostFilterCoeff(single_out_ch_com_sp->data(), single_out_ch_com_sp->data(), post_filter_coeff.data(),
                         kNumFreqs);
}

AMY_AUDIO_API void ComputePowerDensity(float *power_density_prev, float *power_density,
                                       Complexf *ch1_comp_spectrum, Complexf *ch2_comp_spectrum, float smooth_factor,
                                       int len) {
    int i = 0;
    Complexf conj_com_vec = Complexf(0.0f, 0.0f);
    for (i = 0; i < len; i++) {
        conj_com_vec = Complexf(ch2_comp_spectrum[i].real(), -ch2_comp_spectrum[i].imag());
        power_density[i] = smooth_factor * power_density_prev[i] +
                           (1 - smooth_factor) * (ch1_comp_spectrum[i] * conj_com_vec).real();
        power_density_prev[i] = power_density[i];
    }
}


AMY_AUDIO_API void ComputeAveragePowerDensity(float *src_power_density, float *dst_power_density, int chs, int freqs) {
    int i = 0, j = 0, src_chs_ind = 0;
    for (i = 0; i < freqs; i++) {
        for (j = 0; j < chs; j++) {
            src_chs_ind = j * freqs + i;
            dst_power_density[i] += src_power_density[src_chs_ind];
        }
        dst_power_density[i] /= chs;
    }
}

AMY_AUDIO_API void ComputePostFilterCoeff(float *voice_sp, float *noisy_voice_sp, float *filter_coeff, int len) {
    for (int i = 0; i < len; i++) {
        filter_coeff[i] = voice_sp[i] / noisy_voice_sp[i];
    }
}

AMY_AUDIO_API void GetCwiseMaxValOfVecAndScalar(float *vec1, float scalar, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] > scalar ? vec1[i] : scalar;
    }
}

AMY_AUDIO_API void GetCwiseMinValOf2Vecs(float *vec1, float *vec2, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] < vec2[i] ? vec1[i] : vec2[i];
    }
}

AMY_AUDIO_API void GetCwiseMaxValOf2Vecs(float *vec1, float *vec2, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] > vec2[i] ? vec1[i] : vec2[i];
    }
}

AMY_AUDIO_API void GetCwiseMinValOfVecAndScalar(float *vec1, float scalar, float *dst, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = vec1[i] < scalar ? vec1[i] : scalar;
    }
}


AMY_AUDIO_API void
PostFilterGetFilterCoeffMin(float *min_filter_coeff, int low_ind, int high_ind, int len, float slope) {
    int i = 0;
    for (i = 0; i < low_ind - 1; i++)
        min_filter_coeff[i] = 0;
    for (i = low_ind - 1; i < high_ind; i++)
        min_filter_coeff[i] = slope * (i - low_ind + 1) / (high_ind - low_ind);
    for (i = high_ind; i < len; i++)
        min_filter_coeff[i] = slope;

}


AMY_AUDIO_API void
ComputeAverageComSpectrum(Complexf *src_com_spectrum, Complexf *dst_com_spectrum, int chs, int freqs) {
    int i = 0, j = 0, src_chs_ind = 0;
    for (i = 0; i < freqs; i++) {
        dst_com_spectrum[i] = Complexf(0.0f, 0.0f);
        for (j = 0; j < chs; j++) {
            src_chs_ind = j * freqs + i;
            dst_com_spectrum[i] += src_com_spectrum[src_chs_ind];
        }
        dst_com_spectrum[i] /= chs * 1.0f;
    }
}

AMY_AUDIO_API void ConstructComSpFromEsp_Psp(Complexf *dst_sp, float *energy_sp, float *phase_sp, int len) {
    float real_part = 0.0f, imag_part = 0.0f;
    for (int i = 0; i < len; i++) {
        real_part = sqrtf(energy_sp[i]) * cosf(phase_sp[i]);
        imag_part = sqrtf(energy_sp[i]) * sinf(phase_sp[i]);
        dst_sp[i] = Complexf(real_part, imag_part);
    }
}


AMY_AUDIO_API void ApplyPostFilterCoeff(Complexf *src_sp, Complexf *dst_sp, float *filter_coeff, int len) {
    for (int i = 0; i < len; i++) {
        dst_sp[i] = src_sp[i] * filter_coeff[i];
    }
}


AMY_AUDIO_API void FreeMicArrayProc(MicArrayProcInst *self) {
    delete self->mic_array_info;
    self->mic_array_info = NULL;

    delete self->snr_est->smooth_win_of_sp;
    delete self->snr_est->smoothed_local_win;
    delete self->snr_est->smoothed_global_win;
    delete self->snr_est;
    self->snr_est = NULL;


    delete self->vadInst->smooth_win_of_sp;
    delete self->vadInst;
    self->vadInst = NULL;

    delete self->doa_core;
    self->doa_core = NULL;

    delete self->bf_core;
    self->bf_core = NULL;

    delete self->pf_core->min_postfilter_coeff_cn;
    delete self->pf_core;
    self->pf_core = NULL;

}

