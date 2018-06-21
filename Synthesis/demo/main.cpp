#include <iostream>
#include "mic_array_processing.h"
#include "wav_file_header_module.h"
#include "wav_file_module.h"
#include "array_info.h"
#include "audio_file_process.h"
#include <sys/time.h>


using std::cout;
using std::endl;

namespace main_test_single {
    int main(int argc, char *argv[]) {
        FILE *in_fp, *out_fp;

        const int multi_chs_frame_len = kNumChannels * kFrameLen;
        /* input buffer */
        int16_t in_frame_buf_int16[kNumChannels * kFrameLen] = {0};
        float in_frame_buf_float[kNumChannels * kFrameLen] = {0.0f};

        float **in_multi_chs_buf = new float *[kNumChannels];
        for (int i = 0; i < kNumChannels; i++)
            in_multi_chs_buf[i] = new float[kFrameLen]{0.0};

        /* output buffer */
        float out_frame_buf_float[kFrameLen] = {0.0f};
        int16_t out_frame_buf_int16[kFrameLen] = {0};

        /* variable for timing */
        struct timeval start_val, end_val;
        float elapsed_time;

        in_fp = fopen(argv[2], "rb");
        out_fp = fopen(argv[3], "wb");
        if ((NULL == in_fp) || (NULL == out_fp)) {
            cout << "fopen error!" << endl;
            return -1;
        } else {
            cout << "open input file " << argv[1] << endl;
            cout << "open output file " << argv[2] << endl;
        }

        /* variable for file reading and writing operations */
        int read_size = 0, write_size = 0;

        /* instance of MicArrayProcInst */
        MicArrayProcInst pProc;
        const char *array_filename = "array_position_info_2.txt";

        /* init mic array processing instance */
        InitMicArrayProc(&pProc);
        pProc.mic_array_info->azimuth = 0.0f;
        InitMicArrayFromTxtFile(array_filename, pProc.mic_array_info);

        /* set parameters */
        int win_type = 1;
        SetParasMicArrayProc(&pProc, win_type);

        /* validate parameters */
        ValidateParasMicArrayProc(&pProc);

        gettimeofday(&start_val, NULL);
        while (!feof(in_fp)) {
            /* read audio data from input file and do some pre-treatment */
            read_size = fread(in_frame_buf_int16, sizeof(int16_t), multi_chs_frame_len, in_fp);    
            Int16_tToFloat(in_frame_buf_int16, multi_chs_frame_len, in_frame_buf_float);
            Deinterleave(in_frame_buf_float, in_multi_chs_buf);
            /* perform array processing algorithms */
            ProcCoreMicArrayProcSingleOut(&pProc, in_multi_chs_buf, out_frame_buf_float);
            /* write processed audio data into output file */
            FloatToInt16_t(out_frame_buf_float, kFrameLen, out_frame_buf_int16);
            write_size = fwrite(out_frame_buf_int16, sizeof(int16_t), kFrameLen, out_fp);
        }

        gettimeofday(&end_val, NULL);
        elapsed_time = (int64_t)((end_val.tv_sec - start_val.tv_sec) * 1000) +
                       1.f * (end_val.tv_usec - start_val.tv_usec) / 1000.f;

        cout << "Total wasting time is : " << elapsed_time << "(ms)" << endl;
        cout << "Wasting time of each frame is : " << elapsed_time / pProc.frame_num << "(ms)" << endl;

        /* free mic array processing instance */
        FreeMicArrayProc(&pProc);

        fclose(in_fp);
        fclose(out_fp);

        for (int i = 0; i < kNumChannels; i++)
            delete[] in_multi_chs_buf[i];
        delete[] in_multi_chs_buf;

        return 0;
    }
}


namespace main_test_multi {
    int main(int argc, char *argv[]) {
        FILE *in_fp, *out_fp;

        const int multi_chs_frame_len = kNumChannels * kFrameLen;
        /* input buffer */
        int16_t in_frame_buf_int16[multi_chs_frame_len] = {0};
        float in_frame_buf_float[multi_chs_frame_len] = {0.0f};

        float **in_multi_chs_buf = new float *[kNumChannels];
        for (int i = 0; i < kNumChannels; i++)
            in_multi_chs_buf[i] = new float[kFrameLen]{0.0};

        /* output buffer */
        float **out_multi_chs_buf = new float *[kNumChannels];
        for (int i = 0; i < kNumChannels; i++)
            out_multi_chs_buf[i] = new float[kFrameLen]{0};
        float out_frame_buf_float[multi_chs_frame_len] = {0.0f};
        int16_t out_frame_buf_int16[multi_chs_frame_len] = {0};

        /* variable for timing */
        struct timeval start_val, end_val;
        float elapsed_time;

        in_fp = fopen(argv[2], "rb");
        out_fp = fopen(argv[3], "wb");
        if ((NULL == in_fp) || (NULL == out_fp)) {
            cout << "fopen error!" << endl;
            return -1;
        } else {
            cout << "open input file " << argv[1] << endl;
            cout << "open output file " << argv[2] << endl;
        }

        /* variable for file read and write operations */
        int read_size = 0, write_size = 0;

        /* instance of MicArrayProcInst */
        MicArrayProcInst pProc;
        const char *array_filename = "array_position_info_2.txt";

        /* init mic array processing instance */
        InitMicArrayProc(&pProc);
        pProc.mic_array_info->azimuth = 0.0f;
        InitMicArrayFromTxtFile(array_filename, pProc.mic_array_info);

        /* set parameters */
        int win_type = 1;
        SetParasMicArrayProc(&pProc, win_type);

        /* validate parameters */
        ValidateParasMicArrayProc(&pProc);

        gettimeofday(&start_val, NULL);
        while (!feof(in_fp)) {
            /* read audio data from input file and do some pre-treatment */
            read_size = fread(in_frame_buf_int16, sizeof(int16_t), multi_chs_frame_len, in_fp);
            Int16_tToFloat(in_frame_buf_int16, multi_chs_frame_len, in_frame_buf_float);
            Deinterleave(in_frame_buf_float, in_multi_chs_buf);
            /* perform array processing algorithms */
            ProcCoreMicArrayProcMultiOut(&pProc, in_multi_chs_buf, out_multi_chs_buf);
            /* write processed audio data into output file */
            Interleave(out_multi_chs_buf, out_frame_buf_float);
            FloatToInt16_t(out_frame_buf_float, multi_chs_frame_len, out_frame_buf_int16);
            write_size = fwrite(out_frame_buf_int16, sizeof(int16_t), multi_chs_frame_len, out_fp);
        }

        gettimeofday(&end_val, NULL);
        elapsed_time = (int64_t)((end_val.tv_sec - start_val.tv_sec) * 1000) +
                       1.f * (end_val.tv_usec - start_val.tv_usec) / 1000.f;

        cout << "Total wasting time is : " << elapsed_time << "(ms)" << endl;
        cout << "Wasting time of each frame is : " << elapsed_time / pProc.frame_num << "(ms)" << endl;

        /* free mic array processing instance */
        FreeMicArrayProc(&pProc);

        fclose(in_fp);
        fclose(out_fp);

        for (int i = 0; i < kNumChannels; i++)
            delete[] in_multi_chs_buf[i];
        delete[] in_multi_chs_buf;

        for (int i = 0; i < kNumChannels; i++)
            delete[] out_multi_chs_buf[i];
        delete[] out_multi_chs_buf;

        return 0;
    }
}


namespace main_test_single_wav {
    int main(int argc, char *argv[]) {

        const int multi_chs_frame_len = kNumChannels * kFrameLen;

        /* input buffer */
        float in_frame_buf_float16[kNumChannels * kFrameLen] = {0.0f};
        float in_frame_buf_float[kNumChannels * kFrameLen] = {0.0f};

        float **in_multi_chs_buf = new float *[kNumChannels];
        for (int i = 0; i < kNumChannels; i++)
            in_multi_chs_buf[i] = new float[kFrameLen]{0.0};

        /* output buffer */
        float out_frame_buf_float[kFrameLen] = {0.0f};
        float out_frame_buf_float16[kFrameLen] = {0.0f};

        /* variable for timing */
        struct timeval start_val, end_val;
        float elapsed_time;

        const std::string input_file_name(argv[2]);
        const std::string output_file_name(argv[3]);

        WavReader inputfile(input_file_name);
        WavWriter outputfile(output_file_name, inputfile.sample_rate(), 1);

        std::string informat(inputfile.FormatAsString());
        std::string outformat(outputfile.FormatAsString());

        std::cout << "input file format:" << informat << std::endl;
        std::cout << "output file format:" << outformat << std::endl;

        /* instance of MicArrayProcInst */
        MicArrayProcInst pProc;
//	const char *array_filename = "array_position_info_4_cir.txt";

        /* init mic array processing instance */
        InitMicArrayProc(&pProc);
        pProc.mic_array_info->azimuth = 0.0f;
//	InitMicArrayFromTxtFile(array_filename, pProc.mic_array_info);
        InitMicArrayFromFloatArray(XmosArrayPos, pProc.mic_array_info);

        /* set parameters */
        int win_type = 1;
        SetParasMicArrayProc(&pProc, win_type);

        /* validate parameters */
        ValidateParasMicArrayProc(&pProc);


        /* variable for file read and write operations */
        gettimeofday(&start_val, NULL);

        while (inputfile.ReadSamples(multi_chs_frame_len, in_frame_buf_float16) == multi_chs_frame_len) {
            FloatS16ToFloat(in_frame_buf_float16, multi_chs_frame_len, in_frame_buf_float);
            Deinterleave(in_frame_buf_float, in_multi_chs_buf);
            /* perform array processing algorithms */
            ProcCoreMicArrayProcSingleOut(&pProc, in_multi_chs_buf, out_frame_buf_float);
            /* write processed audio data into output file */
            FloatToFloatS16(out_frame_buf_float, kFrameLen, out_frame_buf_float16);
            outputfile.WriteSamples(out_frame_buf_float16, kFrameLen);
        }


        gettimeofday(&end_val, NULL);
        elapsed_time = (int64_t)((end_val.tv_sec - start_val.tv_sec) * 1000) +
                       1.f * (end_val.tv_usec - start_val.tv_usec) / 1000.f;

        cout << "Total wasting time is : " << elapsed_time << "(ms)" << endl;
        cout << "Wasting time of each frame is : " << elapsed_time / pProc.frame_num << "(ms)" << endl;

        /* free mic array processing instance */
        FreeMicArrayProc(&pProc);


        for (int i = 0; i < kNumChannels; i++)
            delete[] in_multi_chs_buf[i];
        delete[] in_multi_chs_buf;

        return 0;
    }
}


int main(int argc, char *argv[]) {

    if (argc < 3) {
        cout << "wrong input parameters: " << argc << endl;
        cout << "Usage: exec_name,exec_path_number,input_file_name,output_file_name" << endl;
        return -1;
    }

    int exec_path_number = atoi(argv[1]);
    if (exec_path_number < 0 || exec_path_number > 2) {
        cout << "value of exec_path_number must be 0,1,or 2" << endl;
        return -1;
    }

    switch (exec_path_number) {
        case 0:
            return main_test_multi::main(argc, argv);
            break;
        case 1:
            return main_test_single::main(argc, argv);
            break;
        case 2:
            return main_test_single_wav::main(argc, argv);
            break;
        default:
            return -1;
    }


}