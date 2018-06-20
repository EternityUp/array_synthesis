#include "wav_file_header_module.h"
#include "wav_file_module.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
    using std::cout;
    using std::endl;
   
    int frame_num = 0, frame_len = 160;
    size_t ret = 0, shift_ret = 0;
    const std::string input_file_name(argv[1]);
    const std::string output_file_name(argv[2]);
    WavReader inputfile(input_file_name);
    WavWriter outputfile(output_file_name, inputfile.sample_rate(), inputfile.num_channels());
    std::string informat(inputfile.FormatAsString());
    std::string outformat(outputfile.FormatAsString());
    cout << "input file format:" << informat << endl;
    cout << "output file format:" << outformat << endl;
    
    size_t unit_size = frame_len * inputfile.num_channels();
    std::vector<int16_t> interleaved(unit_size);
    while ((ret = inputfile.ReadSamples(unit_size, &interleaved[0])) != 0) {
      if (ret < unit_size)
      {
	cout << "end of file" << endl;
	outputfile.WriteSamples(&interleaved[0], ret);
      }
      else
      {
	outputfile.WriteSamples(&interleaved[0], unit_size);
      }
      frame_num++;
      shift_ret += ret;
    }
	
    cout << "frame_num: " << frame_num << "  shift_ret: " << shift_ret << endl; 
    outformat = outputfile.FormatAsString();
    cout << "finally, output file format:" << outformat << endl;
    
    return 0;
}
