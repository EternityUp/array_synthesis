#include "array_process_wrapper.h"

// 构造函数
ArrayProcessSynthesis::ArrayProcessSynthesis() { }

void ArrayProcessSynthesis::InitArrayProcessSynthesis()
{
    printf("Init a array process handle!\n");
    in_chs_audio_data = new float *[kNumChannels];
    for (int i = 0; i < kNumChannels; i++)
      in_chs_audio_data[i] = new float[kFrameLen]{0.0};
    out_ch_audio_data = new float[kFrameLen]{0.0};
    InitMicArrayProc(&process_inst);
}

// 初始化麦克风阵列
void ArrayProcessSynthesis::InitMicArray()
{
    printf("Init a structure containing information of mic array!\n");
    InitMicArrayFromFloatArray(XmosArrayPos, process_inst.mic_array_info);
}


// 设置相关参数
void ArrayProcessSynthesis::SetParas()
{
    printf("Set key parameters of core algorithms!\n");
    int win_type = 1;
    SetParasMicArrayProc(&process_inst, win_type);
}


// 验证参数正确性
void ArrayProcessSynthesis::ValidateParas()
{
    printf("Validate correctness of some key parameters of core algorithms!\n");
    ValidateParasMicArrayProc(&process_inst);
}


// 核心阵列信号处理参数
void ArrayProcessSynthesis::ArrayProcessCore(float *in_multi_chs_data, float *out_single_ch_data)
{
    printf("core of array singal processing algorithms\n");


    for (int i = 0; i < kNumChannels; i++)
    {
      for (int j = 0; j < kFrameLen; j++)
        {
          int index = i * kFrameLen + j;
          in_chs_audio_data[i][j] = in_multi_chs_data[index];
        }
      // printf("in_chs_audio_data[%d][0]=%f\n",i, in_chs_audio_data[i][1]);
    }


    ProcCoreMicArrayProcSingleOut(&process_inst, in_chs_audio_data, out_ch_audio_data);



    for (int i = 0; i < kFrameLen; i++)
    {
      out_ch_audio_data_array[i] = out_ch_audio_data[i];
      // printf("out_ch_audio_data_array[%d]=%f\n", i, out_ch_audio_data[i]);
    }


}


// 释放阵列信号处理算法执行过程中申请的动态内存
void ArrayProcessSynthesis::FreeProcInst()
{
    printf("Free memory allocated dynamically!\n");
    for (int i = 0; i < kNumChannels; i++)
      delete[] in_chs_audio_data[i];
    delete[] in_chs_audio_data;
    delete[] out_ch_audio_data;
    FreeMicArrayProc(&process_inst);
}



