#pragma once 
#include "mic_array_processing.h"


class ArrayProcessSynthesis
{
public:
  ArrayProcessSynthesis();
  void InitArrayProcessSynthesis();
  void InitMicArray();
  void SetParas();
  void ValidateParas();
  void ArrayProcessCore(float **in_multi_chs_data, float *out_single_ch_data);
  void FreeProcInst();
private:
  MicArrayProcInst process_inst;
};




