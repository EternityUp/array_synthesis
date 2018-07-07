#include "array_process_wrapper.h"

// 构造函数
ArrayProcessSynthesis::ArrayProcessSynthesis() { }

void ArrayProcessSynthesis::InitArrayProcessSynthesis()
{
    printf("Init a array process handle!\n");
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
void ArrayProcessSynthesis::ArrayProcessCore(float **in_multi_chs_data, float *out_single_ch_data)
{
    printf("core of array singal processing algorithms\n");
    ProcCoreMicArrayProcSingleOut(&process_inst, in_multi_chs_data, out_single_ch_data);
}


// 释放阵列信号处理算法执行过程中申请的动态内存
void ArrayProcessSynthesis::FreeProcInst()
{
    printf("Free memory allocated dynamically!\n");
    FreeMicArrayProc(&process_inst);
}



