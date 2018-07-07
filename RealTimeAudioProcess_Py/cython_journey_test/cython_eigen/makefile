CPPFLAGS = -g -std=c++11 -I /home/xzc/eigen3.3.4 -I ./inc/ -Wall -fPIC -c

LDFLAGS = -shared -L ./libs/ -lArrayShapeAnalysis -lAudioFileProcess -lFourierTransform -lBasicSignalProcess -lWavFileRW -lm -lstdc++

all:
	g++ $(CPPFLAGS) mic_array_processing.cpp -o mic_array_processing.o
	ld $(LDFLAGS) mic_array_processing.o -o array_process.so


.PHONY:clean
clean:
	rm -f mic_array_processing.o array_process.so
