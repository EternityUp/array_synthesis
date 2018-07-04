# -*- encoding=utf-8 -*-

from ctypes import *

ColVectorXf = POINTER(c_float)

print ColVectorXf(), ColVectorXf


class VaNEstInst(Structure):
    _fields_ = [("smoothed_win_len_of_sp", c_float),
                ("alpha_sp", c_float),
                ("noise_minsp_est_win_len", c_short),
                ("likelihood_ratio_threshold", c_float),
                ("alpha_speech_presence_prob", c_float),
                ("beta_prior_snr", c_float),
                ("smoothed_local_win_len", c_float),
                ("smoothed_global_win_len", c_float),
                ("compromising", c_float),
                ("sigma_min", c_float),
                ("sigma_max", c_float),
                ("sigma_pmin", c_float),
                ("sigma_pmax", c_float),
                ("max_speech_absence_prob", c_float),
                ("Gmin", c_float),
                ("prior_speech_presence", c_float),
                ("kKbdAlpha", c_float),
                ("alpha_ns_prob_update", c_float),
                ("snr_prev_val", c_float),
                ("smooth_win_of_sp", ColVectorXf),
                ("smoothed_local_win", ColVectorXf),
                ("smoothed_global_win", ColVectorXf)]


class VoiceDetInst(Structure):
    _fields_ = [("frequency_start", c_short),
                ("frequency_end", c_short),
                ("frequency_start_ind", c_short),
                ("frequency_end_ind", c_short),
                ("beyond_ratio_energy", c_float),
                ("beyond_ratio_ee", c_float),
                ("initial_voice_prob", c_float),
                ("voice_prob_threshold", c_float),
                ("smoothing", c_float),
                ("smoothed_win_len_of_sp", c_short),
                ("smooth_win_of_sp", ColVectorXf),
                ("LargerNum", c_int),
                ("alpha_pre_emphasis", c_float),
                ("alpha_log_sp", c_float),
                ("speech_start_threshold", c_float),
                ("speech_finished_threshold", c_float),
                ("non_speech_to_speech_frame_len", c_ushort),
                ("speech_to_non_speech_frame_len", c_ushort),
                ("is_speech_frame", c_bool),
                ("is_speech_triggered", c_bool)]


kNumFreqs = 129
freq_in_array = (c_float * kNumFreqs)


class DoaInst(Structure):
    _fields_ = [("frequency_start", c_float),
                ("frequency_end", c_float),
                ("frequency_start_ind", c_short),
                ("frequency_end_ind", c_short),
                ("freq_bin", freq_in_array),
                ("coarse_ang_span", c_float),
                ("coarse_ang_start", c_float),
                ("coarse_ang_max1", c_float),
                ("coarse_ang_max2", c_float),
                ("coarse_beam_num", c_float),
                ("fine_center_ang", c_float),
                ("fine_ang_span", c_float),
                ("loading_factor", c_float),
                ("sp_smoothing", c_float),
                ("doa", c_float)]


class bfInst(Structure):
    _fields_ = [("freq_bin", freq_in_array)]


class PfInst(Structure):
    _fields_ = [("Num_Inters", c_int),
                ("smooth_alpha", c_float),
                ("min_postfilter_coeff", c_float),
                ("freq_low", c_float),
                ("freq_high", c_float),
                ("freq_low_ind", c_int),
                ("freq_high_ind", c_int),
                ("slope", c_float),
                ("min_postfilter_coeff_cn", ColVectorXf)]


class Point(Structure):
    _fields_ = [("", c_float),
                ("", c_float),
                ("", c_float)]


class ArrayProps(Structure):
    _fields_ = [("location",c_int),
                ("array_normal", ),
                ("array_direction",),
                ("min_mic_spacing",),
                ("mic_spacings",),
                ("is_linear",),
                ("is_planar",),
                ("is_tridimensional",),
                ("element_num",),
                ("azimuth",),
                ("sound_path",)]


class MicArrayProcInst(Structure):
    _fields_ = [("mic_array_info", ),
                ("snr_est",),
                ("vadInst", ),
                ("vad_mode", c_int),
                ("doa_core", ),
                ("doa_mode", c_int),
                ("bf_core", ),
                ("bf_mode", c_int),
                ("pf_core", ),
                ("cn_mode", c_int),
                ("buf_rw_status", ),
                ("block_tf", ),
                ("window", ),
                ("noise_sps_est",),
                ("recur_avg_sps_est",),
                ("prior_snr",),
                ("post_snr",),
                ("gain",),
                ("speech_prob", ),
                ("is_speech", c_bool),
                ("is_speech_triggered", c_bool),
                ("doa", c_float),
                ("frame_num", c_uint),
                ("frame_offset", c_int),
                ("kKbdAlpha", c_float),
                ("freq_resolution", c_float),]

